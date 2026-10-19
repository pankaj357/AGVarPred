/** @type {import('tailwindcss').Config} */
export default {
  content: ["./index.html", "./src/**/*.{ts,tsx}"],
  theme: {
    extend: {
      colors: {
        ink: "#0b2020",
        panel: "#102a29",
        mint: "#a2e6c6",
        aqua: "#61d7cc",
      },
      boxShadow: {
        glow: "0 0 0 1px rgba(162,230,198,.18), 0 12px 40px rgba(0,0,0,.2)",
      },
    },
  },
  plugins: [],
};
