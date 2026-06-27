document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('.summary-note, .model-subnote, .summary-table').forEach(function (el) {
    el.remove();
  });
  document.querySelectorAll('.summary-pills .tf-tag + .tf-tag').forEach(function (el) {
    el.remove();
  });
  function cleanModelText(text) {
    text = text.split(' · ')[0];
    text = text.split('active PDBs').join('3D models');
    text = text.split('active PDB models').join('3D models');
    text = text.split('active PDB').join('3D model');
    if (text.indexOf('3D model') === -1) {
      text = text.split(' models').join(' 3D models');
      text = text.split(' model').join(' 3D model');
    }
    text = text.split('3D 3D models').join('3D models');
    text = text.split('3D 3D model').join('3D model');
    return text;
  }
  document.querySelectorAll('.model-region summary span:last-child').forEach(function (el) {
    el.textContent = cleanModelText(el.textContent);
  });
  document.querySelectorAll('.region-card dd').forEach(function (el) {
    var text = cleanModelText(el.textContent);
    text = text.split('Predicted = Low').join('Predicted');
    el.textContent = text;
  });
  var walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  var nodes = [];
  while (walker.nextNode()) nodes.push(walker.currentNode);
  nodes.forEach(function (node) {
    var text = node.nodeValue;
    text = text.split('Predicted = Low').join('Predicted');
    text = text.split('Collapsible technical table. Rows are grouped by protein region and retain model, template, residue coverage, and available summary evidence.').join('');
    text = text.split('active PDB models').join('3D models');
    text = text.split('active PDB model').join('3D model');
    text = text.split('3D 3D models').join('3D models');
    text = text.split('3D 3D model').join('3D model');
    text = text.split('Scan this region').join('Scan');
    node.nodeValue = text;
  });
});
