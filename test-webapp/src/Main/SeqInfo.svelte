<script lang="ts">
  import { Peptidoglycan } from "smithereens";
  export let value;
</script>

<div class="flex gap-3 font-plex-mono">
  <span><strong>Monoisotopic Mass: </strong>{value.monoisotopic_mass()}</span>
</div>
