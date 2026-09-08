import test from 'node:test';
import assert from 'node:assert/strict';
import {coordinate,extent,format,ticks} from '../live/plot.mjs';
test('logarithmic coordinates keep small positive objective gaps',()=>{
  assert.equal(coordinate(1e-12,[1e-15,1e-3],0,120,true),30);
  assert.deepEqual(ticks([1e-5,1e-2],true),[1e-5,1e-4,1e-3,1e-2]);
});
test('constant and empty series still have usable boxed axes',()=>{
  assert.ok(extent([2,2])[1]>extent([2,2])[0]);
  assert.ok(extent([-2,-2])[1]>extent([-2,-2])[0]);
  assert.deepEqual(extent([]),[0,1]);
  assert.ok(extent([0,-1,NaN,1e-12],true)[0]>0);
  assert.ok(extent([1e-20,1e-18],true)[0]<1e-20);
});
test('linear ticks use readable interior values without changing the data bounds',()=>{
  assert.deepEqual(ticks([-1.2,21.2]),[0,5,10,15,20]);
  assert.deepEqual(ticks([-.06,1.06]),[0,.25,.5,.75,1]);
});
test('readouts do not round tiny nonzero values to zero',()=>{
  assert.equal(format(1e-15),'1.0e-15');assert.equal(format(0),'0');
});
