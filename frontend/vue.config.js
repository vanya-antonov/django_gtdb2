module.exports = {
  transpileDependencies: ["vuetify"],
  publicPath: "/",
  configureWebpack: {
    module: {
      rules: [
        {
          test: /\.js$/,
          use: [
            "ify-loader",
            "transform-loader?plotly.js/tasks/compress_attributes.js",
          ],
        },
      ],
    },
  },
};
