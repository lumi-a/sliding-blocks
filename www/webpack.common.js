const CopyWebpackPlugin = require("copy-webpack-plugin");
const path = require("path");

module.exports = {
  entry: "./bootstrap.ts",
  output: {
    path: path.resolve(__dirname, "dist"),
    filename: "bootstrap.js",
    clean: true
  },
  resolve: {
    extensions: [".ts", ".js", ".wasm"],
  },
  module: {
    rules: [
      {
        test: /\.ts$/,
        use: "ts-loader",
        exclude: /node_modules/
      }
    ]
  },
  plugins: [new CopyWebpackPlugin(["index.html", "style.css", "favicon.ico"])],
  experiments: {
    asyncWebAssembly: true,
    syncWebAssembly: true,
  },
};