/* Run with Playwright and Chromium. Set CHROMIUM_PATH for a system browser. */
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { pathToFileURL } = require("node:url");
const { chromium } = require("playwright");
const base = path.resolve(__dirname, "..");
const dist = fs.existsSync(path.join(base, "dist")) ? path.join(base, "dist") : base;

(async () => {
  const browser = await chromium.launch({
    executablePath: process.env.CHROMIUM_PATH || undefined, args: ["--no-sandbox"]
  });
  try {
    const page = await browser.newPage();
    const errors = [];
    page.on("pageerror", e => errors.push(e.message));
    await page.goto(pathToFileURL(path.join(dist, "structure.html")).href);
    const run = async (button) => {
      await page.locator(button).click();
      await page.waitForFunction(() => document.querySelector("#status").textContent !== "Yapı işleniyor…");
    };
    const record = await page.evaluate(() => window.polymerReference);
    await run("#reference");
    assert.equal(await page.locator("#geometry").isVisible(), true);
    assert.equal(await page.locator('#next a[href="contact.html"]').count(), 1);

    await page.locator('[data-tab="import"]').click();
    const setFile = async obj => {
      await page.locator("#import-file").setInputFiles({ name:"structure.json", mimeType:"application/json", buffer:Buffer.from(JSON.stringify(obj)) });
      await page.waitForFunction(() => document.querySelector("#status").textContent !== "Yapı işleniyor…");
    };
    const bad = structuredClone(record);
    bad.structure.canonical_smiles = "CC";
    await setFile(bad);
    assert.match(await page.locator("#status").innerText(), /eşleşmiyor/);
    assert.equal(await page.locator("#output").isVisible(), false);
    await setFile(record);
    assert.equal(await page.locator("#geometry").isVisible(), true);

    await page.locator('[data-tab="smiles"]').click();
    const input = async smiles => {await page.locator("#smiles").fill(smiles);await run("#build");};
    await input("C(C)(C)(C)(C)C");
    assert.match(await page.locator("#status").innerText(), /değerliği aşıyor/);
    await input("[13CH3][C@H](O)C");
    assert.equal(await page.locator("#geometry").isVisible(), true);
    assert.equal(await page.locator('#next a[href="contact.html"]').count(), 0);
    const href = await page.locator("#next a").first().getAttribute("href");
    assert.match(href, /system=unsupported/);
    const downloadEvent = page.waitForEvent("download");
    await page.locator("#save").click();
    const download = await downloadEvent;
    const saved = JSON.parse(fs.readFileSync(await download.path(), "utf8"));
    assert.equal(saved.input_smiles, "[13CH3][C@H](O)C");
    await page.locator('[data-tab="import"]').click();
    await setFile(saved);
    assert.equal(await page.locator("#geometry").isVisible(), true);
    assert.equal(href, await page.locator("#next a").first().getAttribute("href"));

    await page.locator('[data-tab="smiles"]').click();
    await input("CC=CC");
    assert.equal(await page.locator("#geometry").isVisible(), false);
    await input(record.structure.canonical_smiles);
    assert.equal(await page.locator('#next a[href="contact.html"]').count(), 1);

    await page.locator('[data-tab="sketch"]').click();
    await page.locator('#sketch [data-action="add"]').click();
    await page.locator('#sketch [data-action="add"]').click();
    await run("#build-sketch");
    assert.match(await page.locator("#identity").innerText(), /C2H6/);
    assert.equal(await page.locator("#geometry").isVisible(), true);
    assert.deepEqual(errors, []);
    console.log("Structure import, identity checks, stereochemistry, atom editor, 3D geometry and export round-trip passed.");
  } finally { await browser.close(); }
})().catch(err => {console.error(err);process.exitCode=1;});
