// Interpretation outlet only: all atmospheric, geometric and spectral
// coefficients originate in the Python exporter. No ray solver lives here.
export async function readJSON(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
  return response.json();
}

async function digest(bytes) {
  const hash = await crypto.subtle.digest('SHA-256', bytes);
  return Array.from(new Uint8Array(hash), b => b.toString(16).padStart(2, '0')).join('');
}

export function validateManifest(m) {
  if (m.schema_version !== 1) throw new Error('Unsupported asset schema');
  if (!m.common || !m.families || !m.presets?.opac || !m.domain) throw new Error('Incomplete manifest');
  if (m.data_bytes > 32*2**20) throw new Error('Asset budget exceeded');
  for (const item of [m.common, ...Object.values(m.families)]) {
    if (!/^data\/[a-z0-9-]+\.bin\.gz$/.test(item.file) || item.encoding !== 'gzip' ||
        !/^[a-f0-9]{64}$/.test(item.sha256) || !Number.isSafeInteger(item.decodedByteLength) ||
        item.decodedByteLength <= 0 || item.decodedByteLength > 64*2**20) throw new Error('Invalid asset descriptor');
    const ranges=[];
    for (const field of Object.values(item.fields)) {
      const width=field.dtype === '<u4' ? 4 : field.dtype === '<f8' ? 8 : 0;
      const count=field.shape.reduce((a,b)=>a*b,1);
      if (!width || field.shape.some(n=>!Number.isSafeInteger(n)||n<0) ||
          !Number.isSafeInteger(count) || field.byteLength !== count*width ||
          !Number.isSafeInteger(field.offset) || field.offset<0 || field.offset%width ||
          field.offset+field.byteLength>item.decodedByteLength) throw new Error('Invalid array shape/offset');
      ranges.push([field.offset,field.offset+field.byteLength]);
    }
    ranges.sort((a,b)=>a[0]-b[0]);
    if(ranges.some((x,i)=>i && x[0]<ranges[i-1][1])) throw new Error('Overlapping arrays');
  }
  for(const preset of Object.values(m.presets)) {
    if(!m.families[preset.family] || preset.amplitudes.length!==preset.layers.length+2 ||
       preset.amplitudes.some(x=>!Number.isFinite(x)||x<0)) throw new Error('Invalid preset/layer pairing');
  }
  return m;
}

async function fetchPacked(base, descriptor) {
  const response=await fetch(new URL(descriptor.file,base));
  if(!response.ok) throw new Error(`${descriptor.file}: HTTP ${response.status}`);
  const buffer=await response.arrayBuffer();
  if(buffer.byteLength!==descriptor.byteLength || await digest(buffer)!==descriptor.sha256) throw new Error(`${descriptor.file}: asset checksum mismatch`);
  return buffer;
}

export async function decode(buffer, descriptor) {
  // The bytes are explicitly gzip, not dependent on HTTP Content-Encoding.
  const raw=await new Response(new Blob([buffer]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
  if(raw.byteLength!==descriptor.decodedByteLength) throw new Error('Decoded asset size mismatch');
  if(new Uint8Array(new Uint32Array([1]).buffer)[0]!==1) throw new Error('Little-endian platform required');
  const arrays={};
  for(const [name,f] of Object.entries(descriptor.fields)) {
    const Type=f.dtype==='<u4'?Uint32Array:Float64Array;
    arrays[name]=new Type(raw,f.offset,f.byteLength/Type.BYTES_PER_ELEMENT);
    for(const value of arrays[name]) if(!Number.isFinite(value)) throw new Error(`${name}: nonfinite asset`);
  }
  return arrays;
}

export async function loadSite(base=new URL('./',import.meta.url)) {
  const manifest=validateManifest(await readJSON(new URL('manifest.json',base)));
  const common=await decode(await fetchPacked(base,manifest.common),manifest.common);
  if(common.wave.length!==81 || common.radius[0]!==0 || common.radius.at(-1)!==9000000 ||
     common.cie.length!==243 || common.rgbMatrix.length!==9 || common.reference.length!==81 ||
     common.referenceXYZ.length!==3 || common.referenceXYZ[1]<=0 ||
     common.radius.some((v,i)=>i && v<=common.radius[i-1])) throw new Error('Invalid common axes');
  const data={manifest,common,families:{}};
  // Four byte-deduplicated families fit the explicit memory budget. Decode
  // sequentially once; discard compressed buffers instead of repeatedly paying
  // decompression/validation on every interactive preset switch.
  for(const [key,descriptor] of Object.entries(manifest.families)) {
    data.families[key]=await decode(await fetchPacked(base,descriptor),descriptor);
    validateFamily(data,key,data.families[key]);
  }
  data.workspace=new Float64Array(Math.max(...Object.values(data.families).map(a=>a.b.length))*common.wave.length);
  return data;
}

function validateFamily(data,key,a) {
  const descriptor=data.manifest.families[key];
  const nr=data.common.radius.length,nw=data.common.wave.length,nc=descriptor.fields.components.shape[1];
  const nn=a.b.length;
  if(a.offsets.length!==nr+1 || a.outer.length!==nr || a.index.length!==a.weight.length ||
     a.offsets[0]!==0 || a.offsets.at(-1)!==a.index.length ||
     a.components.length!==nn*nc || a.factors.length!==nc*nw ||
     a.offsets.some((v,i)=>i && v<a.offsets[i-1]) || a.index.some(v=>v>=nn) ||
     a.weight.some(v=>v<0) || a.components.some(v=>v<0) || a.outer.some(v=>v<0)) throw new Error('Invalid receiver/component pairing');
  for(const p of Object.values(data.manifest.presets)) {
    if(p.family===key && p.amplitudes.length!==nc) throw new Error('Preset component count mismatch');
  }
}

export async function selectFamily(data,key) {
  const arrays=data.families[key];
  if(!arrays)throw new Error('Unknown family');
  return arrays;
}

export function checkParameters(m, parameters) {
  const {preset,du,aod,rh}=parameters;
  if(!Object.hasOwn(m.presets,preset)) throw new Error('Unknown preset');
  for(const [name,value] of Object.entries({du,aod,rh})) {
    const [lo,hi]=m.domain[name];
    if(!Number.isFinite(value)||value<lo||value>hi) throw new Error(`${name}: outside supported UI domain`);
  }
  if(aod>0 && aod<m.domain.positive_aod_min) throw new Error('AOD below minimum positive step');
}

export async function reconstruct(data, parameters) {
  checkParameters(data.manifest,parameters);
  const {manifest:m,common:c}=data;
  const selected=m.presets[parameters.preset];
  // CLI zero override disables the *whole* column, including ST. It is not
  // applied to fixed calibrated presets or none, which do not allow overrides.
  const scalable=parameters.preset==='opac'||parameters.preset==='opac-ft';
  const preset=scalable && parameters.aod===0?m.presets.none:selected;
  const a=await selectFamily(data,preset.family);
  const amplitudes=preset.amplitudes.slice();
  amplitudes[1]=parameters.du;
  if(scalable && parameters.aod>0) {
    const knots=m.rh_knots,i=parameters.rh<=knots[1][0]?0:1;
    const growth=knots[i][1]+(parameters.rh-knots[i][0])/(knots[i+1][0]-knots[i][0])*(knots[i+1][1]-knots[i][1]);
    amplitudes[2]=parameters.aod*growth;
  }
  const nw=c.wave.length,nr=c.radius.length,nc=amplitudes.length,nn=a.b.length;
  const nodes=data.workspace;
  for(let i=0;i<nn;i++) {
    for(let k=0;k<nw;k++) {
      let tau=a.components[i*nc]*a.factors[k]+a.components[i*nc+1]*parameters.du*a.factors[nw+k];
      for(let j=2;j<nc;j++) tau+=a.components[i*nc+j]*amplitudes[j]*a.factors[j*nw+k];
      nodes[i*nw+k]=Math.exp(-tau);
    }
  }
  const {offsets,index:ids,weight:weights,outer}=a;
  const transmission=new Float64Array(nr*nw);
  for(let r=0;r<nr;r++) {
    const dest=r*nw,lo=offsets[r],hi=offsets[r+1];
    let k=0;
    // Scalar blocking shares exact index/weight loads. It never merges nodes,
    // drops annuli, interpolates exp(tau), or changes spectral integration order.
    for(;k+7<nw;k+=8) {
      let s0=0,s1=0,s2=0,s3=0,s4=0,s5=0,s6=0,s7=0;
      for(let j=lo;j<hi;j++) {
        const src=ids[j]*nw+k,w=weights[j];
        s0+=w*nodes[src];s1+=w*nodes[src+1];s2+=w*nodes[src+2];s3+=w*nodes[src+3];
        s4+=w*nodes[src+4];s5+=w*nodes[src+5];s6+=w*nodes[src+6];s7+=w*nodes[src+7];
      }
      transmission[dest+k]=s0+outer[r];transmission[dest+k+1]=s1+outer[r];
      transmission[dest+k+2]=s2+outer[r];transmission[dest+k+3]=s3+outer[r];
      transmission[dest+k+4]=s4+outer[r];transmission[dest+k+5]=s5+outer[r];
      transmission[dest+k+6]=s6+outer[r];transmission[dest+k+7]=s7+outer[r];
    }
    for(;k<nw;k++) {
      let sum=0;
      for(let j=lo;j<hi;j++) sum+=weights[j]*nodes[ids[j]*nw+k];
      transmission[dest+k]=sum+outer[r];
    }
  }
  const xyz=new Float64Array(nr*3),rgb=new Float64Array(nr*3),relativeY=new Float64Array(nr);
  for(let r=0;r<nr;r++) {
    for(let k=0;k<nw;k++) {
      const e=transmission[r*nw+k]*c.reference[k];
      for(let j=0;j<3;j++) xyz[r*3+j]+=e*c.cie[k*3+j];
    }
    for(let j=0;j<3;j++) {
      for(let q=0;q<3;q++) rgb[r*3+j]+=xyz[r*3+q]*c.rgbMatrix[j*3+q];
    }
    for(let j=0;j<3;j++) xyz[r*3+j]/=c.referenceXYZ[1];
    relativeY[r]=xyz[r*3+1];
  }
  return {transmission,xyz,rgb,relativeY};
}

export function interval(radius, r) {
  if(!Number.isFinite(r)||r<radius[0]||r>radius.at(-1)) throw new Error('Radius outside exported domain');
  let lo=0,hi=radius.length-1;
  while(hi-lo>1) {const mid=(lo+hi)>>1;if(radius[mid]<=r)lo=mid;else hi=mid;}
  return [lo,(r-radius[lo])/(radius[lo+1]-radius[lo])];
}

export function sample(values, radius, r, channels=3) {
  const [i,w]=interval(radius,r);
  return Array.from({length:channels},(_,c)=>values[i*channels+c]*(1-w)+values[(i+1)*channels+c]*w);
}

export function finish(rgb,ev) {
  if(!Number.isFinite(ev)||ev<-4||ev>20) throw new Error('Exposure outside UI domain');
  const gain=2**ev,values=rgb.map(x=>Math.max(0,x)*gain),scale=Math.max(1,...values);
  return values.map(x=>{const v=x/scale;return v<=.0031308?12.92*v:1.055*v**(1/2.4)-.055;});
}

export function pixelMap(radius,moonRadius,x,y,size=560,previous=null) {
  if(!Number.isFinite(x)||!Number.isFinite(y)||Math.hypot(x,y)>radius.at(-1)-moonRadius+1e-6) throw new Error('Moon center outside supported domain');
  const index=previous?.size===size?previous.index:new Int32Array(size*size);
  const weight=previous?.size===size?previous.weight:new Float64Array(size*size);
  index.fill(-1);
  for(let py=0;py<size;py++) for(let px=0;px<size;px++) {
    const dx=((px+.5)/size*2-1)*moonRadius,dy=((py+.5)/size*2-1)*moonRadius;
    if(dx*dx+dy*dy>moonRadius*moonRadius)continue;
    const r=Math.hypot(dx+x,dy+y);
    let lo=0,hi=radius.length-1;
    while(hi-lo>1){const mid=(lo+hi)>>1;if(radius[mid]<=r)lo=mid;else hi=mid;}
    index[py*size+px]=lo;weight[py*size+px]=(r-radius[lo])/(radius[lo+1]-radius[lo]);
  }
  return {index,weight,size,image:previous?.size===size?previous.image:null};
}

export function draw(ctx,rgb,map,ev) {
  const size=map.size,image=map.image??ctx.createImageData(size,size),bytes=image.data,gain=2**ev;
  map.image=image;
  const encode=v=>255*(v<=.0031308?12.92*v:1.055*v**(1/2.4)-.055);
  for(let p=0;p<map.index.length;p++) {
    const i=map.index[p]*3;if(i<0)continue;
    const w=map.weight[p],u=1-w;
    const red=Math.max(0,rgb[i]*u+rgb[i+3]*w)*gain;
    const green=Math.max(0,rgb[i+1]*u+rgb[i+4]*w)*gain;
    const blue=Math.max(0,rgb[i+2]*u+rgb[i+5]*w)*gain;
    const scale=Math.max(1,red,green,blue),offset=p*4;
    bytes[offset]=encode(red/scale);bytes[offset+1]=encode(green/scale);bytes[offset+2]=encode(blue/scale);bytes[offset+3]=255;
  }
  ctx.putImageData(image,0,0);
  return image;
}
