#!/usr/bin/env node
// scrnaseq-lab: the Lab's converter, from a shell or a pipeline.
//
// The app's own modules run by Node instead of a browser tab: each part is the
// bundle the app's builder writes, and what the command adds (sentinels, its
// notes, a sha256 per part) is listed in the README. The app is unchanged by
// this: it still converts in the browser and uploads nothing.
// See src/cli/main.ts, or run `scrnaseq-lab --help`.
//
// This file is plain JavaScript on purpose. The modules it runs are TypeScript,
// which Node runs directly by stripping the types (Node 22.18 and later) — but
// only when this build of Node can, and only outside node_modules. Both are
// checked here first, so that either failure is a sentence and not a trace.

if (!process.features?.typescript) {
  process.stderr.write(`scrnaseq-lab: Node ${process.version} cannot run TypeScript sources; this command needs Node 22.18 or later\n`)
  process.exit(1)
}

let main
try {
  ({ main } = await import('../src/cli/main.ts'))
} catch (e) {
  if (e?.code === 'ERR_UNSUPPORTED_NODE_MODULES_TYPE_STRIPPING') {
    process.stderr.write('scrnaseq-lab: this copy is inside node_modules, where Node will not run TypeScript. Use a clone: git clone https://github.com/JiaenLin/scrnaseq-lab && cd scrnaseq-lab && npm ci && npm link\n')
    process.exit(1)
  }
  throw e
}
process.exitCode = await main(process.argv.slice(2))
