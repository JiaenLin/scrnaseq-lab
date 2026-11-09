/** @type {import('tailwindcss').Config} */
export default {
  // Everything that renders. Not src/cli: the command draws nothing, and a word
  // in its comments ('grow', 'block', 'hidden') would otherwise ship as CSS.
  content: ['./index.html', './src/*.{ts,tsx}', './src/components/**/*.{ts,tsx}', './src/lib/**/*.{ts,tsx}'],
  theme: {
    extend: {
      fontFamily: {
        mono: ['ui-monospace', 'SFMono-Regular', 'Consolas', 'monospace'],
      },
    },
  },
  plugins: [],
}
