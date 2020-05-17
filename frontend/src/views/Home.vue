<template>
    <v-data-table
        :headers="headers"
        :items="organisms"
        item-key="name"
        fixed-header
        disable-pagination
        show-group-by
        hide-default-footer
        height="calc(100vh - 130px)"
        style="max-height: calc(100vh - 144px)"
    >
        <template v-slot:item.name="{ item }">
    <router-link
        :to="{ name: 'OrganismDetails', params: { id: item.id } }"
        target="_blank"
        >{{ item.name }}</router-link
    >
</template>
    </v-data-table>
</template>
<script>
import HTTP from "@/http-common";

export default {
    data: function() {
        return {
            headers: [
                { text: "Organism", value: "name" },
                { text: "Phylum", value: "phylum" },
                { text: "Kingdom", value: "kingdom" },
                { text: "Genotype", value: "genotype" },
            ],
            organisms: [],
            height: 100,
        };
    },
   

    created: function() {
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
