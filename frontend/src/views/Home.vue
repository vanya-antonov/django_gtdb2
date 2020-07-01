<style>
/*#nohideticks .y2tick2 {
	display: none !important;
}
#nohideticks .y2tick {
	display: block !important;*/
/*}/*/
</style>

<template>
		<div >
			<div class="mx-auto" style="position: right;">
			<h1>Chelatase DB</h1>
			<v-card max-width="30em">
				Magnesium chelatase chlIDH and cobalt chelatase cobNST enzymes are required for biosynthesis of (bacterio)chlorophyll and cobalamin (vitamin B12), respectively. Each enzyme consists of large, medium and small subunits. Structural and primary sequence similarities indicate common evolutionary origin of the corresponding subunits. This also allows prokaryotes to use taxa-specific sets of genes to encode the magnesium and/or cobalt chelatases.
In our recent study we analyzed 1200+ prokaryotic genomes (see the Organisms section) and identified at least 5 taxa-specific strategies that can be used to synthesize these two chelatases (Antonov 2020):
Interestingly, more than 100 species may utilize programmed frameshifting in the chlD gene to produce both the small and the medium chelatase subunits from the same gene (the predicted frameshifting signals can be found in the Signals section)
Below is a graphical representation of the 5 proposed strategies.
For visualization purposes, a chelatase is depicted as a sausage (small subunit) and a hot dog (medium subunit) on a plate (large subunit). This representation highlights the facts that each of the chelatases consists of three different subunits and that the small subunit is similar to a part of the medium subunit (in other words, “a sausage is a part of a hot dog”).

			</v-card>
		</div>
			<chelTable ></chelTable>
			<plotlyHeatmap
				v-if="organisms.length"
				:click-callback="clickCallback"
				:organisms="organisms"
				:height="700"
			></plotlyHeatmap>
			<br /><br/>
			<v-row justify="right"><v-spacing></v-spacing>
			<h2 class="mx-auto" style="position: right;">
				<router-link
					v-if="clickedIndex"
					:to="{
						name: 'OrganismDetails',
						params: { id: organisms[clickedIndex].id },
					}"
					target="_blank"
					>{{ organisms[clickedIndex].name }}</router-link
				>
			</h2>
</v-row ><v-spacing></v-spacing>
			<!-- <v-row id="nohideticks" class="d-flex justify-center"> -->
				<plotlyHeatmap
					v-if="clickedIndex"
					:organisms="zoomedOrganisms"
					:drawBar="true"
					:click-callback="smallClickCallback"
					:removeY2ticks="false"
					:addLinksToOrganism="true"
				></plotlyHeatmap>
				<!-- </v-row> -->
		</div>
</template>
<script>
import HTTP from "@/http-common";
import PlotlyHeatmap from "@/components/PlotlyHeatmap";
import ChelTable from "@/components/ChelTable";

export default {
	components: { PlotlyHeatmap, ChelTable },
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
			console.log(eventData);
			// if (eventData.points[0].data.type == "heatmap") {
			// 	this.clickedIndex = eventData.points[0].pointIndex[0];
			// }
		},
		setZoomedOrganisms() {
			this.zoomedOrganisms = this.organisms.slice(
				Math.max(0, this.clickedIndex - 3),
				Math.min(this.clickedIndex + 4, this.organisms.length)
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
