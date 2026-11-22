// Add a "Copy" button to each code block.
// Lines that start with "#>" (R output) are not copied, like on pkgdown sites.
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll("div.highlighter-rouge, figure.highlight").forEach(function (block) {
    var code = block.querySelector("pre code") || block.querySelector("pre");
    if (!code) return;

    var button = document.createElement("button");
    button.className = "copy-code-button";
    button.type = "button";
    button.textContent = "Copy";
    button.setAttribute("aria-label", "Copy code to clipboard");

    button.addEventListener("click", function () {
      var text = code.innerText
        .split("\n")
        .filter(function (line) { return !/^#>/.test(line); })
        .join("\n")
        .replace(/\n+$/, "");

      navigator.clipboard.writeText(text).then(
        function () { button.textContent = "Copied"; },
        function () { button.textContent = "Error"; }
      );
      setTimeout(function () { button.textContent = "Copy"; }, 1500);
    });

    block.style.position = "relative";
    block.appendChild(button);
  });
});
