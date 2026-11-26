// Run with Node, Playwright, and Chromium. The test uses the distributable data.
'use strict';
const fs=require('node:fs');
const os=require('node:os');
const path=require('node:path');
const assert=require('node:assert/strict');
const {pathToFileURL}=require('node:url');
const {chromium}=require('playwright');
(async()=>{
  const base=path.resolve(__dirname,'..');
  const web=process.env.TRANSFER_APP_DIR || (fs.existsSync(path.join(base,'dist','index.html'))?path.join(base,'dist'):base);
  const browser=await chromium.launch({executablePath:process.env.CHROMIUM_PATH||undefined,args:['--no-sandbox']});
  const context=await browser.newContext({viewport:{width:390,height:844},acceptDownloads:true});
  try{
    const page=await context.newPage(),errors=[];
    page.on('pageerror',err=>errors.push(String(err)));
    await page.goto(pathToFileURL(path.join(web,'transfer.html')).href);
    await page.waitForSelector('#finding:visible');
    assert.match(await page.locator('#finding').innerText(),/20 eşfaz ölçümünün 4 tanesinde/);
    assert.match(await page.locator('#finding').innerText(),/0,042 ve 0,0455/);
    assert.equal(await page.locator('#plots svg').count(),2);
    assert.equal(await page.locator('#observation option').count(),20);
    await page.locator('#model').selectOption('combined');
    assert.equal(await page.locator('#legend .swatch.cloud').count(),0);
    await page.locator('#observation').selectOption({index:3});
    await page.locator('#model').selectOption('all');
    assert.equal(await page.locator('#scores tbody tr').count(),2);
    await page.locator('[data-view="osmotic"]').click();
    assert.equal(await page.locator('#observation option').count(),3);
    assert.equal(await page.locator('#plots svg').count(),1);
    await page.locator('#sample').selectOption('@calibration');
    assert.equal(await page.locator('#observation option').count(),5);
    await page.locator('#tempUnit').selectOption('K');
    assert.match(await page.locator('#observation option').first().innerText(),/308,15 K/);
    await page.locator('[data-view="cloud"]').click();
    assert.equal(await page.locator('#observation option').count(),24);
    await page.locator('#sample').selectOption('@calibration');
    assert.equal(await page.locator('#observation option').count(),73);
    const dl=page.waitForEvent('download');
    await page.locator('#download').click();
    const file=await dl;assert.equal(file.suggestedFilename(),'PS_siklohekzan_cloud.csv');
    const dir=fs.mkdtempSync(path.join(os.tmpdir(),'transfer-csv-'));
    try {
      const csvPath=path.join(dir,'cloud.csv');await file.saveAs(csvPath);
      const csv=fs.readFileSync(csvPath,'utf8').trim().split(/\r?\n/);
      assert.equal(csv.length,74);
      assert.match(csv[0],/polymer_volume_fraction/);
    } finally {fs.rmSync(dir,{recursive:true,force:true});}
    await page.locator('[data-view="coex"]').click();
    assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth),true);
    assert.deepEqual(errors,[]);
    console.log('PASS · eşfaz 20 çift, iki model, sıcaklık birimi, A₂ dış/uyum 3/5, bulutlanma dış/uyum 24/73, CSV ve mobil genişlik');
  } finally {await browser.close();}
})().catch(err=>{console.error(err);process.exit(1)});
