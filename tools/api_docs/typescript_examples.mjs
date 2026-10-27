#!/usr/bin/env node
// Checks and runs a repository's TypeScript example programs.
//
//  1. Type-checks the examples with the package's own `tsc`.
//  2. Runs each example with Node's type stripping, resolving the package's npm
//     name (e.g. `@earthsciml/ast`) to the in-repo build. Any failure fails the run.
//
// The region markers are checked by regions.py, which typescript.sh runs first.
//
// Usage: node typescript_examples.mjs <repo-root> <pkg-dir> <examples-dir> <npm-name> <dist-entry>
//   pkg-dir and examples-dir are relative to <repo-root>; dist-entry is relative
//   to pkg-dir. typescript.sh reads them from docs/api/config.toml.

import { spawnSync } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'

import { fileURLToPath } from 'node:url'

const repoRoot = path.resolve(process.argv[2] ?? '.')
const [pkgRel, exRel, npmName, distEntry] = process.argv.slice(3)
if (!pkgRel || !exRel || !npmName || !distEntry) {
  console.error('usage: typescript_examples.mjs <repo-root> <pkg-dir> <examples-dir> <npm-name> <dist-entry>')
  process.exit(2)
}
const exDir = path.join(repoRoot, exRel)
const pkgDir = path.join(repoRoot, pkgRel)
const hook = path.join(path.dirname(fileURLToPath(import.meta.url)), 'typescript_examples_resolve.mjs')
const hookEnv = {
  ...process.env,
  API_DOCS_TS_NAME: npmName,
  API_DOCS_TS_ENTRY: path.join(pkgDir, distEntry),
}
const rel = (p) => path.relative(repoRoot, p)

const examples = fs
  .readdirSync(exDir)
  .filter((f) => f.endsWith('.ts') && !f.endsWith('.d.ts'))
  .sort()
  .map((f) => path.join(exDir, f))

let failed = false
const fail = (msg) => {
  console.error(msg)
  failed = true
}

// --- 1. type-check -------------------------------------------------------------
const tsc = path.join(pkgDir, 'node_modules', 'typescript', 'bin', 'tsc')
const tc = spawnSync(process.execPath, [tsc, '-p', path.join(exDir, 'tsconfig.json')], { stdio: 'inherit' })
if (tc.status !== 0) fail('type-check of the TypeScript examples failed')

// --- 2. run --------------------------------------------------------------------
const nodeFlags = ['--disable-warning=ExperimentalWarning']
if (!process.features?.typescript) nodeFlags.push('--experimental-strip-types')
for (const file of examples) {
  const t0 = Date.now()
  const r = spawnSync(process.execPath, [...nodeFlags, '--import', hook, file], {
    stdio: 'inherit',
    cwd: repoRoot,
    env: hookEnv,
  })
  const ms = Date.now() - t0
  if (r.status === 0) console.log(`ok    ${rel(file)} (${ms} ms)`)
  else fail(`FAIL  ${rel(file)} (exit ${r.status ?? r.signal})`)
}

if (failed) {
  console.error('typescript examples: FAILED')
  process.exit(1)
}
console.log(`typescript examples: ${examples.length} passed`)
