/** Every suggested comparison must resolve to a real computed experiment. */
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {getGuide,controlHelp,displayValue} from '../live/guidance.mjs';
const folder=new URL('../live/data/',import.meta.url);
const experiments=fs.readdirSync(folder).filter(name=>name.endsWith('.json')).map(name=>JSON.parse(fs.readFileSync(new URL(name,folder))));
test('all interactive figures have meaningful guidance and executable comparisons',()=>{
  assert.equal(experiments.length,145);
  for(const experiment of experiments){
    const guide=getGuide(experiment);
    assert.ok(guide.prompt.length>35&&guide.observe.length>50,experiment.id);
    assert.ok(guide.presets.length>=2&&guide.presets.length<=3,experiment.id);
    const keys=new Set();
    for(const preset of guide.presets){
      const indices=experiment.controls.map(c=>c.default);
      for(const [id,value] of Object.entries(preset.values)){
        const j=experiment.controls.findIndex(c=>c.id===id);assert.ok(j>=0,`${experiment.id}: ${id}`);
        indices[j]=experiment.controls[j].values.indexOf(value);assert.ok(indices[j]>=0,`${experiment.id}: ${id}=${value}`);
      }
      const key=indices.join(',');assert.ok(experiment.states[key],`${experiment.id}: ${key}`);keys.add(key);
    }
    assert.equal(keys.size,guide.presets.length,`${experiment.id}: duplicate suggested states`);
  }
});
test('parameter explanations distinguish dilation, analysis, and display-only changes',()=>{
  assert.match(controlHelp({view:'wavelets',family:'wavelet1d'},{id:'level'}),/contracts/);
  assert.match(controlHelp({view:'wavcoefs-2d',family:'wavelet2d'},{id:'level'}),/coarser/);
  assert.match(controlHelp({view:'discrete-operators',family:'gradient'},{id:'gain'}),/unchanged/);
  assert.equal(displayValue('db4'),'Daubechies 4');assert.equal(displayValue(0),'0');
});
