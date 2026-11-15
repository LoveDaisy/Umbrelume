import {checkParameters,draw,pixelMap,sample} from './kernel.js';

const $=id=>document.getElementById(id);
const canvas=$('moon'),ctx=canvas.getContext('2d',{willReadFrequently:true});
const worker=new Worker(new URL('./worker.js',import.meta.url),{type:'module'});
const state={preset:'opac',du:300,aod:.013,rh:.5,ev:8,x_m:2860000,y_m:0};
let manifest,common,field,map,active=null,pending=null,sequence=0,latest=0,fatal=false;
let readyResolve,readyReject;
const ready=new Promise((resolve,reject)=>{readyResolve=resolve;readyReject=reject;});
// Initialization errors are displayed, not silently replaced with stale data.
ready.catch(()=>{});
const timing=[];

function fail(message) {
  fatal=true;
  $('status').textContent=`无法计算：${message}。请重新加载；若问题持续，请核对资产与浏览器支持。`;
  $('status').classList.add('error');
  $('controls').disabled=true;
  for(const item of [active,pending])item?.reject(new Error(message));
  active=pending=null;readyReject(new Error(message));worker.terminate();
}

function labels() {
  for(const key of ['du','aod','rh','ev']) {
    const value=state[key];
    $(key).value=key==='aod'?Math.round(value*1000):key==='rh'?Math.round(value*100):value;
    $(`${key}-value`).textContent=key==='du'?`${value} DU`:key==='aod'?value.toFixed(3):key==='rh'?`${Math.round(value*100)}%`:`${value>=0?'+':''}${value.toFixed(1).replace('.0','')} EV`;
  }
  $('preset').value=state.preset;
  const preset=manifest?.presets[state.preset],scalable=state.preset==='opac'||state.preset==='opac-ft';
  $('aod').disabled=!scalable;
  $('rh').disabled=!scalable||state.aod===0;
  $('identity').classList.toggle('fitted',preset?.kind==='calibrated');
  $('identity').textContent=preset?.kind==='calibrated'?'已拟合，非物理先验。固定事件校准态，不接受 AOD / RH 覆盖。':
    state.preset==='none'?'分子大气反事实，无气溶胶。臭氧仍独立可调。':
    state.aod===0?'精确零态：FT 与平流层气溶胶均关闭。':'固定文献先验；湿度只缩放 FT 幅度。平流层背景独立保持 0.005。';
  $('position').textContent=`月心距影轴 ${(Math.hypot(state.x_m,state.y_m)/1000).toFixed(0)} km`;
  $('moon-marker').setAttribute('cx',state.x_m/1000);
  $('moon-marker').setAttribute('cy',state.y_m/1000);
}

function readout() {
  if(!field||!document.querySelector('.readout').open)return;
  const r=Math.hypot(state.x_m,state.y_m),t=sample(field.transmission,common.radius,r,81),y=sample(field.relativeY,common.radius,r,1)[0];
  $('photometry').textContent=`月心相对明视觉照度 Eᵥ / Fᵥ = ${y.toExponential(5)}；数值是径向场的线性插值，不含显示曝光。`;
  $('spectra').replaceChildren(...Array.from(common.wave,(wave,i)=>{
    const tr=document.createElement('tr');
    for(const text of [wave.toFixed(0),(t[i]*common.reference[i]).toExponential(6),t[i].toExponential(6)]) {
      const td=document.createElement('td');td.textContent=text;tr.append(td);
    }
    return tr;
  }));
}

function paint() {
  if(!field)return;
  const start=performance.now();
  if(!map)map=pixelMap(common.radius,manifest.moon_radius_m,state.x_m,state.y_m);
  draw(ctx,field.rgb,map,state.ev);
  readout();
  return performance.now()-start;
}

function dispatch() {
  if(active||!pending||fatal)return;
  active=pending;pending=null;
  worker.postMessage({type:'compute',id:active.id,parameters:active.parameters});
}

function update(parameters={}) {
  if(!manifest||fatal)return Promise.reject(new Error('Demo is not ready'));
  const merged={preset:state.preset,du:state.du,aod:state.aod,rh:state.rh,...parameters};
  checkParameters(manifest,merged);
  Object.assign(state,merged);labels();
  const id=++sequence;latest=id;
  // Bounded queue: one executing request and one latest pending request.
  if(pending)pending.resolve({superseded:true,id:pending.id});
  const promise=new Promise((resolve,reject)=>{pending={id,parameters:merged,start:performance.now(),resolve,reject};});
  dispatch();return promise;
}

function display(values={}) {
  const next={...state,...values},limit=9000000-manifest.moon_radius_m;
  if(!Number.isFinite(next.ev)||next.ev<-4||next.ev>20||!Number.isFinite(next.x_m)||!Number.isFinite(next.y_m)||Math.hypot(next.x_m,next.y_m)>limit+1e-6)throw new Error('Display parameters outside supported domain');
  if(state.x_m!==next.x_m||state.y_m!==next.y_m)map=null;
  state.ev=next.ev;state.x_m=next.x_m;state.y_m=next.y_m;
  labels();return paint();
}

worker.onerror=event=>fail(event.message);
worker.onmessage=({data:message})=>{
  if(message.type==='error'){fail(message.error);return;}
  if(message.type==='ready') {
    manifest=message.manifest;common=message.common;
    $('preset').replaceChildren(...Object.entries(manifest.presets).map(([key,preset])=>{
      const option=document.createElement('option');option.value=key;
      option.textContent=key==='opac'?'OPAC 双层文献先验':key==='opac-ft'?'OPAC 仅自由对流层':key==='none'?'分子大气（无气溶胶）':`${preset.provenance.event}（已拟合）`;
      return option;
    }));
    const [umbra,penumbra]=manifest.boundaries_m.map(x=>x/1000);
    $('umbra').setAttribute('r',umbra);$('penumbra').setAttribute('r',penumbra);
    $('umbra-label').setAttribute('y',-umbra-350);$('penumbra-label').setAttribute('y',-penumbra-350);
    $('version').textContent=`模型 v${manifest.software_version}；资产 schema ${manifest.schema_version}；${(manifest.data_bytes/2**20).toFixed(2)} MiB。首次完整预载后才开放状态切换。`;
    map=pixelMap(common.radius,manifest.moon_radius_m,state.x_m,state.y_m);
    update().then(()=>{$('controls').disabled=false;labels();readyResolve();}).catch(error=>fail(error.message));
    return;
  }
  if(message.type!=='field'||!active||message.id!==active.id)return;
  const request=active;active=null;
  if(message.id===latest) {
    field=message.field;
    const drawTime=paint(),total=performance.now()-request.start;
    const measurement={compute:message.compute,draw:drawTime,total,communication:Math.max(0,total-message.compute-drawTime),parameters:request.parameters};
    timing.push(measurement);if(timing.length>250)timing.shift();
    $('status').textContent='';
    request.resolve({superseded:false,id:message.id,timing:measurement});
  } else request.resolve({superseded:true,id:message.id});
  dispatch();
};

for(const key of ['preset','du','aod','rh'])$(key).addEventListener('input',()=>{
  const value=key==='preset'?$(key).value:Number($(key).value)/(key==='aod'?1000:key==='rh'?100:1);
  update({[key]:value}).catch(error=>fail(error.message));
});
$('ev').addEventListener('input',()=>display({ev:Number($('ev').value)}));
$('reset').addEventListener('click',()=>{
  display({ev:8,x_m:2860000,y_m:0});
  update({preset:'opac',du:300,aod:.013,rh:.5}).catch(error=>fail(error.message));
});
document.querySelector('.readout').addEventListener('toggle',readout);

function move(x,y) {
  if(!manifest||fatal)return;
  const limit=9000000-manifest.moon_radius_m,length=Math.hypot(x,y),scale=length>limit?limit/length:1;
  display({x_m:x*scale,y_m:y*scale});
}
for(const element of [canvas,$('navigator')]) {
  let drag;
  element.addEventListener('pointerdown',event=>{
    if(!manifest||fatal)return;
    event.preventDefault();element.focus();element.setPointerCapture(event.pointerId);
    drag={px:event.clientX,py:event.clientY,x:state.x_m,y:state.y_m};
    if(element!==canvas) {
      const rect=element.getBoundingClientRect();
      move(((event.clientX-rect.left)/rect.width*19000-9500)*1000,((event.clientY-rect.top)/rect.height*19000-9500)*1000);
      drag={px:event.clientX,py:event.clientY,x:state.x_m,y:state.y_m};
    }
  });
  element.addEventListener('pointermove',event=>{
    if(!drag)return;
    const scale=element===canvas?18000000/element.clientWidth:19000000/element.clientWidth;
    move(drag.x+(event.clientX-drag.px)*scale,drag.y+(event.clientY-drag.py)*scale);
  });
  for(const type of ['pointerup','pointercancel','lostpointercapture'])element.addEventListener(type,()=>{drag=null;});
  element.addEventListener('keydown',event=>{
    const direction={ArrowLeft:[-1,0],ArrowRight:[1,0],ArrowUp:[0,-1],ArrowDown:[0,1]}[event.key];
    if(!direction||!manifest)return;
    event.preventDefault();const step=event.shiftKey?500000:50000;
    move(state.x_m+direction[0]*step,state.y_m+direction[1]*step);
  });
}

// The validation harness exercises this same application/worker/canvas path;
// these are inspectable results, not an alternative test-only implementation.
window.demo={ready,update,display,get state(){return {...state};},get field(){return field;},
  get manifest(){return manifest;},get common(){return common;},get timing(){return timing.slice();},
  get queue(){return {active:active?.id??null,pending:pending?.id??null,latest};}};
worker.postMessage({type:'initialize'});
