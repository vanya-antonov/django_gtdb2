<style>
.ytick {
	display: none !important;
	width: 0;
}
</style>

<template>
	<v-row>
		<v-col>
			<plotlyHeatmap style="height: 1000; width: 500;"></plotlyHeatmap
		></v-col>
	</v-row>
</template>
<script>
import PlotlyHeatmap from "@/components/PlotlyHeatmap";
export default {
	components: { PlotlyHeatmap },
};
</script>
