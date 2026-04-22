window.MathJax = {
  loader: {
    load: ['[tex]/autoload', '[tex]/upgreek', '[tex]/unicode', '[tex]/cancel', '[tex]/boldsymbol', '[tex]/newcommand', '[tex]/color'],
  },
  tex: {
    packages: {'[+]': ['autoload', 'upgreek', 'unicode', 'cancel', 'boldsymbol', 'newcommand', 'color']},
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$'], ['\\[', '\\]']],
    processEscapes: true,
    processEnvironments: true,
    macros: {
      RR: '{\\mathbb{R}}',
      bold: ['{\\bf #1}', 1],
      eo: '{\\varepsilon_0}',
      e: ['{\\times 10^{#1}}', 1],
      epsilon: '{\\varepsilon}',
      k: '{1.381\\times10^{-23}}',
      q: '{1.602\\times10^{-19}}',
      eon: '{8.854\\times10^{-12}}',
      cnmr: '{\\ce{^13C-NMR}}',
      hnmr: '{\\ce{^1H-NMR}}',
      AA: '{\\unicode{x212B}}',
      infin: '{\\unicode{x221E}}',
      kjmol: '{KJ\\cdot mol^{-1}}',
      il: ['{[\\ce{#1}][\\ce{#2}]}', 2]
    }
  },
  options: {
    ignoreHtmlClass: '.*|',
    processHtmlClass: 'arithmatex'
  }
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
