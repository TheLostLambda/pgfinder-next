<script lang="ts">
  import { Peptidoglycan, pg_to_fragments } from "smithereens";
  import { writable } from "svelte/store";
  import Grid from "./assets/Grid.svg";
  import Input from "./Main/Input.svelte";
  import Output from "./Main/Output.svelte";

  let input = "gm-AEJA";
  let pipeline = writable((s) => s.to_string());
  let output = "";
  let error = false;

  let swap = () => (input = output);
  let copy = () => navigator.clipboard.writeText(output);

  let seq = new Peptidoglycan(input);
  $: try {
    // FIXME: Move whitespace handling to the Rust side of things!
    seq = new Peptidoglycan(input.replace(/\s/g, ""));
    error = false;
  } catch (e) {
    output = e;
    error = true;
  }

  $: try {
    if (!error) {
      output = pg_to_fragments(seq);
    }
  } catch (e) {
    output = e;
  }
</script>

<main
  class="flex portrait:flex-col landscape:flex-row justify-evenly items-center bg-surface-200 dark:bg-surface-800 bg-repeat"
  style="background-image: url({Grid})"
>
  <Input bind:value={input} bind:seq />
  <Output bind:value={output} on:swap={swap} on:copy={copy} />
</main>
