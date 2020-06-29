<style>
#nohideticks .y2tick2 {
	display: none !important;
}
#nohideticks .y2tick {
	display: block !important;
}
</style>

<template>
	<v-col cols="auto">
		<chelTable></chelTable>
		<div fluid class="d-flex justify-end">
			<plotlyHeatmap
				v-if="organisms.length"
				:click-callback="clickCallback"
				:organisms="organisms"
				:height="700"
			></plotlyHeatmap>
		</div>
		<router-link
			v-if="clickedIndex"
			:to="{
				name: 'OrganismDetails',
				params: { id: organisms[clickedIndex].id },
			}"
			target="_blank"
			>{{ organisms[clickedIndex].name }}</router-link
		>
		<v-row id="nohideticks" class="d-flex justify-end">
			<plotlyHeatmap
				v-if="clickedIndex"
				:organisms="zoomedOrganisms"
				:drawBar="true"
				:click-callback="smallClickCallback"
			></plotlyHeatmap>
		</v-row>
		<router-link :to="{name: 'OrganismsList',
		params: {kingdomsToFilterProp: ['Bacteria'],
		phylumsToFilterProp: ['Actinobacteria', 'Proteobacteria']}}">Organisms 222</router-link>
	</v-col>
</template>
<script>
import HTTP from "@/http-common";
import PlotlyHeatmap from "@/components/PlotlyHeatmap";
import ChelTable from "@/components/ChelTable";

export default {
	components: { PlotlyHeatmap,
		ChelTable},
	data: function() {
		return { organisms: [], clickedIndex: null, zoomedOrganisms: [] };
	},
	methods: {
		clickCallback(eventData) {
			console.log("callback");
			console.log(eventData);
			if (eventData.points[0].data.type == "heatmap") {
				this.clickedIndex = eventData.points[0].pointIndex[0];
				this.setZoomedOrganisms();
			}
		},
		smallClickCallback(eventData) {
			console.log(eventData)
			// if (eventData.points[0].data.type == "heatmap") {
			// 	this.clickedIndex = eventData.points[0].pointIndex[0];
			// }
		},
		setZoomedOrganisms() {
			this.zoomedOrganisms = this.organisms.slice(
				Math.max(0, this.clickedIndex - 10),
				Math.min(this.clickedIndex + 10, this.organisms.length)
			);
			console.log(this.zoomedOrganisms);
		},
	},

	created() {
		HTTP.get(`organisms`)
			.then((response) => {
				// JSON responses are automatically parsed.
				this.organisms = response.data;
			})
			.catch((error) => {
				console.warn(error);
			});
	},
};
</script>
