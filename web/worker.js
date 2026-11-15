import {loadSite,reconstruct} from './kernel.js';

let data;
self.onmessage=async({data:message})=>{
  try {
    if(message.type==='initialize') {
      data=await loadSite(new URL('./',import.meta.url));
      self.postMessage({type:'ready',manifest:data.manifest,common:data.common});
      return;
    }
    if(message.type!=='compute'||!data)throw new Error('Worker is not initialized');
    const start=performance.now();
    const field=await reconstruct(data,message.parameters);
    const compute=performance.now()-start;
    self.postMessage({type:'field',id:message.id,parameters:message.parameters,field,compute},
      Object.values(field).map(a=>a.buffer));
  } catch(error) {
    self.postMessage({type:'error',id:message.id,error:error.message});
  }
};
