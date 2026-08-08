import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// base './' → works under GitHub Pages /<repo>/ with no extra config.
export default defineConfig({
  base: './',
  plugins: [react()],
  // The reader worker imports h5wasm on demand, which is a code split, and
  // Vite's default IIFE worker format cannot express one.
  worker: { format: 'es' },
})
