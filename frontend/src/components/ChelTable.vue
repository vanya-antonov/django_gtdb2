<style>
#cheltable th {
	font-size: 1.5rem !important;
}
#cheltable td {
	max-width: 500px !important;
	white-space: nowrap;
}
</style>
<template
	>
		<v-simple-table id="cheltable">
			<template v-slot:default>
				<thead>
					<th class="text-h4" style="width: 400px;">
						Production of Mg-/Co-chelatases
					</th>
					<th class="text-h4">Possible Genotype(s)</th>
				</thead>
				<tbody>
					<tr key="0">
						<td>
							<v-img
								position="left"
								:height="130"
								contain
								:src="require('@/assets/images/chelTable0.png')"
							/>
						</td>
						<td>
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xchlH',
													'1xchlD',
													'1xchlI',
												],
												type: 'exact',
											},
											phylumsToFilterProp: [
												'Proteobacteria',
											],
										}),
									},
								}"
								target="_blank"
								><p class="font-weight-bold">
									Proteobacteria: chlH, chlD, chl
								</p></router-link
							>
							<p>Mg-chelatase: (Gibson et al. 1995)</p>
						</td>
					</tr>
					<tr key="1">
						<td>
							<v-img
								position="left"
								:height="130"
								contain
								:src="require('@/assets/images/chelTable1.png')"
							/>
						</td>
						<td>
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xcobN',
													'1xcobT',
													'1xcobS',
													'1xchlH',
													'1xchlD',
													'1xchlI',
												],
												type: 'exact',
											},
											phylumsToFilterProp: [
												'Proteobacteria',
											],
										}),
									},
								}"
								target="_blank"
								><p class="font-weight-bold">
									Proteobacteria: cobN, cobT, cobS, chlH,
									chlD, chlI
								</p></router-link
							>
							<p>
								Mg-chelatase: (Gibson et al. 1995)<br />
								Co-chelatase: (Debussche et al. 1992)
							</p>
						</td>
					</tr>
					<tr key="2">
						<td>
							<v-img
								position="left"
								:height="130"
								contain
								:src="require('@/assets/images/chelTable2.png')"
							/>
						</td>
						<td>
							<p class="font-weight-bold">
								<router-link
									:to="{
										name: 'OrganismsList',
										query: {
											filters: JSON.stringify({
												genesToFilterProp: {
													values: [
														'1xcobN',
														'1xchlH',
														'1xchlD',
														'1xchlI',
													],
													type: 'exact',
												},
												phylumsToFilterProp: [
													'Cyanobacteria',
												],
											}),
										},
									}"
									target="_blank"
								>
									Cyanobacteria: cobN, chlH, chlD, chlI </router-link
								><br />
								<router-link
									:to="{
										name: 'OrganismsList',
										query: {
											filters: JSON.stringify({
												genesToFilterProp: {
													values: [
														'1xcobN',
														'2xchlH',
														'1xchlD',
														'1xchlI',
													],
													type: 'exact',
												},
												phylumsToFilterProp: [
													'Cyanobacteria',
												],
											}),
										},
									}"
									target="_blank"
								>
									Cyanobacteria: cobN, 2xchlH, chlD, chlI
								</router-link>
							</p>
							<p>
								Mg-chelatase: (Gibson et al. 1995)<br />Co-chelatase:
								(Rodionov et al. 2003)
							</p>
						</td>
					</tr>
					<tr key="3">
						<td>
							<v-img
								position="left"
								:height="130"
								contain
								:src="require('@/assets/images/chelTable3.png')"
							/>
						</td>
						<td>
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xcobN',
													'1xchlD',
													'1xchlI',
												],
												type: 'exact',
											},
											phylumsToFilterProp: [
												'Actinobacteria',
											],
										}),
									},
								}"
								target="_blank"
								><p class="font-weight-bold">
									Actinobacteria: cobN, chlD, chlI
								</p></router-link
							>
							<p>Co-chelatase: (Rodionov  et al. 2003)</p>
						</td>
					</tr>
					<tr key="4">
						<td>
							<v-img
								position="left"
								:height="130"
								contain
								:src="require('@/assets/images/chelTable4.png')"
							/>
						</td>
						<td><p class="font-weight-bold">
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xcobN', '1xfs-chlD',
												],
												type: 'exact',
											},
											phylumsToFilterProp: [
												'Proteobacteria',
											],
										}),
									},
								}"
								target="_blank"
								>
									Proteobacteria: cobN, fs-chlD
								</router-link
							><br/>
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xcobN', '1xfs-chlD',
												],
												type: 'exact',
											},
											phylumsToFilterProp: [
												'Actinobacteria',
											],
										}),
									},
								}"
								target="_blank"
								>
									Actinobacteria: cobN, fs-chlD
								</router-link
							><br/>
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xcobN', '1xfs-chlD',
												],
												type: 'exact',
											},
											kingdomsToFilterProp: [
												'Archaea',
											],
										}),
									},
								}"
								target="_blank"
								>
									Archaea: cobN, fs-chlD
								</router-link
							>
							</p>
							<p>Co-chelatase: (Antonov 2020)</p>
						</td>
					</tr>
				<!-- 	<tr key="5">
						<td>
							<v-img
								position="left"
								:height="130"
								contain
								:src="require('@/assets/images/chelTable5.png')"
							/>
						</td>
						<td>
							<p class="font-weight-bold">Bold text.</p>
							<router-link
								:to="{
									name: 'OrganismsList',
									query: {
										filters: JSON.stringify({
											genesToFilterProp: {
												values: [
													'1xchlH',
													'1xchlD',
													'1xchlI',
												],
												type: 'exact',
											},
											phylumsToFilterProp: [
												'Proteobacteria',
											],
										}),
									},
								}"
								target="_blank"
								><p class="font-weight-bold">
									Proteobacteria: chlH, chlD, chlI
								</p></router-link
							>
							<p>reference</p>
						</td>
					</tr> -->
				</tbody>
			</template>
		</v-simple-table>
</template>
<script>
export default {};
</script>
