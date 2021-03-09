<template>
    <v-app>
        <v-app-bar absolute color="white" elevation="0" app>
            <v-btn text large :to="'/'" retain-focus-on-click="true">
                <v-toolbar-title>Chelatase DB</v-toolbar-title>
            </v-btn>
            <v-btn text large :to="'/organisms'">Organisms</v-btn>
            <v-btn text large :to="'/signals'">Signals</v-btn>
            <v-btn text large :to="'/about'">About</v-btn>
        </v-app-bar>
        <v-content> 
            <v-container fluid>
                <router-view />
            </v-container>
        </v-content>
    </v-app>
</template>
<script>
// import HelloWorld from './components/HelloWorld';

export default {
    name: "App",
    data: () => ({
        //
    }),
};
</script>
