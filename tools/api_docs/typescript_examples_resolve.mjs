// Node `--import` hook for the TypeScript examples: resolves the package's npm
// name (API_DOCS_TS_NAME, e.g. `@earthsciml/ast`) to the in-repo build
// (API_DOCS_TS_ENTRY, an absolute path), so the examples import the package
// exactly as a user would but run against this checkout rather than a
// published release. typescript_examples.mjs sets both variables.
import module from 'node:module'
import { pathToFileURL } from 'node:url'

const entry = pathToFileURL(process.env.API_DOCS_TS_ENTRY).href
const NAME = process.env.API_DOCS_TS_NAME

if (typeof module.registerHooks === 'function') {
  module.registerHooks({
    resolve: (specifier, context, next) =>
      specifier === NAME ? { url: entry, shortCircuit: true } : next(specifier, context),
  })
} else {
  const src =
    `export async function resolve(s, c, n) {` +
    ` return s === ${JSON.stringify(NAME)} ? { url: ${JSON.stringify(entry)}, shortCircuit: true } : n(s, c) }`
  module.register('data:text/javascript,' + encodeURIComponent(src))
}
