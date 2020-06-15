<template>
    <div ref="fornaContainer" class="d-inline-flex pa-2"></div>
</template>
<script>
import { FornaContainer as PreFornaContainer } from "fornac";

export default {
    name: "fornac",
    data: function() {
        return {
            fornaContainer: [],
            options: {
                structure: "((..((....)).(((....))).))",
                sequence: "CGCUUCAUAUAAUCCUAAUGACCUAU",
                labelInterval: 10,
                transitionDuration: 500,
                friction: 30,
            },
            colorScheme: "structure",
        };
    },

    mounted: function() {
        this.fornaContainer[0] = new PreFornaContainer(this.$refs.fornaContainer, {
            initialSize: [600, 600],
            allowPanningAndZooming: true,
        });
        this.fornaContainer[0].changeColorScheme(this.colorScheme);
        this.fornaContainer[0].addRNA(this.options.structure, this.options);
        this.fornaContainer[0].center_view();
    },
};
</script>
