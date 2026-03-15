// MathJax v3 configuration for MkDocs Material
// Place this file BEFORE the MathJax CDN script in mkdocs.yml

window.MathJax = {
  tex: {
    inlineMath: [['\\(', '\\)'], ['$', '$']],    // inline math
    displayMath: [['\\[', '\\]'], ['$$', '$$']], // display math
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    processHtmlClass: 'arithmatex', // only process .arithmatex spans
    ignoreHtmlClass: '.*'
  },
  startup: {
    typeset: false // disable auto-typeset; we'll call it ourselves
  }
};

// Hook into MkDocs Material instant navigation
if (window.document$) {
  document$.subscribe(() => {
    const container = document.querySelector('.md-content');
    if (window.MathJax) {
      // Ensure MathJax has finished startup before typesetting
      MathJax.startup.promise
        .then(() => MathJax.typesetPromise([container]))
        .catch(err =>
          console.error("MathJax typeset failed: " + err.message)
        );
    }
  });
}
