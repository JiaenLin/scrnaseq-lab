import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// base './' → works under GitHub Pages /<repo>/ with no extra config.
export default defineConfig({
  base: './',
  plugins: [react()],
  // The reader is a module worker; h5wasm itself is loaded at runtime from
  // public/vendor so the bundler never rewrites it.
  worker: { format: 'es' },
})
