<style>
#fstext {
    min-width: 23rem;
}
#organismsabout {
    min-width: 35rem;
}
</style>
<template>
    <div>
    <v-row justify="center" class=" px-6">
        <v-col id="fstext">
            <h1 style=" white-space: nowrap;">Programed Frameshifting</h1>
            <p>
            <ul>
                <li>
                    Translational or programmed ribosomal frameshifting (PRF) is
                    a special type of recoding utilized for expression of some
                    genes and can be found in all domains of life. On PRF a
                    ribosome changes the initial reading frame at a specific
                    location in mRNA. Our previous large-scale computational
                    study of more than 1100 prokaryotic genomes has identified a
                    number of putative genes with PRF (<a href="https://pubmed.ncbi.nlm.nih.gov/23649834/" target="_blank">Antonov et al. 2013</a>).
                    Subsequent experimental verification has confirmed
                    translational frameshifting in 4 out of 20 gene families
                    with the chelatase chlD gene demonstrating the highest
                    efficiency (up to 60%).
                </li>
                <li>
                    Particularly, in the 2013 study we predicted and
                    experimentally validated highly efficient -1 frameshifting
                    signals in the magnesium chelatase chlD gene (<a href="https://pubmed.ncbi.nlm.nih.gov/23649834/" target="_blank">Antonov et al. 2013</a>):
                </li>
                <ul>
                    <li>
                        The frameshifts in the chlD genes were predicted in <a href="http://83.149.211.146:23194/~ivan/cgi-bin/GeneTack/cgi/cof_view.cgi?id=131733585" target="_blank">23
                        diverse prokaryotes</a>
                    </li>
                    <li>
                        <a href="http://83.149.211.146:23194/~ivan/cgi-bin/GeneTack/cgi/print_page.cgi?fn=conserved_fs.html&title=Conserved%20Frameshifts" target="_blank">Three frameshifting signals</a> were experimentally
                        validated
                    </li>
                </ul>
                <li>
                    Our most recent study (<a href="https://pubmed.ncbi.nlm.nih.gov/32211852/" target="_blank">Antonov 2020</a>) has shown that the PRF
                    mechanism allows a single chlD gene to produce both the
                    small and the medium subunits of the aerobic cobalt
                    chelatase. In that study we analyzed more than 1200
                    prokaryotic genomes and identified more than 150
                    frameshifted chlD genes from both archaea and bacteria.
                    Importantly many of these genes contained putative recoding
                    signals.
                </li>
                <li></li>
                <li>
                    Particularly, we observed several different strategies that
                    can be used by various prokaryotic taxa to produce these
                    enzymes.
                </li>
                <li>
                    The Chelatase DB provides access and allows to browse all
                    the data obtained in the <a href="https://pubmed.ncbi.nlm.nih.gov/32211852/" target="_blank">Antonov 2020</a>.
                </li>
            </ul>
            </p>
        </v-col>
        <v-col>
            <v-img
                position="center"
                :height="700"
                contain
                :src="require('@/assets/images/aboutImage.png')"
            />
        </v-col>
    </v-row>
    <v-spacing/>
    <v-row id="organismsabout" justify="center" class=" px-6">
        <v-spacer/>
        <v-col><h1 align="center">About: <router-link to="/organisms" target="_blank">Organisms Page</router-link></h1>
            <p>A list of all the prokaryotic genomes from the Reference and Representative NCBI sets containing the chlD gene(s). The Mg- and Co-chelatases can be synthesized from several different gene sets (“<b>chelatase genotypes</b>” -- see the last column). Namely, the three component magnesium chelatase can be encoded by a full set of three separate genes (the “<b>chlH, chlD, chlI</b>” genotype). Similarly, the cobalt chelatase cobNST can be produced from the “<b>cobN, cobT, cobS</b>” genotype. Additionally, the small and the medium Mg-chelatase subunits can also function as parts of the chimeric Co-chelatase (the "<b>cobN, chlD, chlI</b>" and "<b>cobN, chlH, chlD, chlI</b>" genotypes). Finally, the translational frameshifting may allow chelatase synthesis from a reduced two-gene genotype (the "<b>cobN, fs-chlD</b>" genotype).
You can use the search form or the filter panel on the right for your convenience.
</p>
        </v-col>
        <v-spacer/>
    </v-row>
</div>
</template>
<script>
    export default{
        mounted() {
            var section=this.$router.currentRoute.hash.replace("#", "");
            if (section)
                this.$nextTick(()=> window.document.getElementById(section).scrollIntoView());}
    }
</script>