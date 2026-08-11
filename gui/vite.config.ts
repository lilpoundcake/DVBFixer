import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'
import { apiPlugin } from './server/api-plugin'

export default defineConfig({
  plugins: [
    react(),
    apiPlugin(),
  ],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
  optimizeDeps: {
    include: [
      'style-to-js',
      'style-to-object',
      'hast-util-to-jsx-runtime',
      'react-markdown',
    ],
  },
  ssr: {
    noExternal: ['molstar'],
  },
})
