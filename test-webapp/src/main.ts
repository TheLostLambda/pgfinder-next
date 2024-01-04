import "@brainandbones/skeleton/themes/theme-vintage.css";
import "@brainandbones/skeleton/styles/all.css";
import "./app.postcss";
import App from "./App.svelte";
import init from "smithereens";

const app = (async () => {
  await init();
  new App({
    target: document.body,
  });
})();

export default app;
