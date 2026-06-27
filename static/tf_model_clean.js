document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('.summary-note, .model-subnote, .summary-table').forEach(function (el) {
    el.style.display = 'none';
  });
  document.querySelectorAll('.summary-pills .tf-tag + .tf-tag').forEach(function (el) {
    el.style.display = 'none';
  });
  document.querySelectorAll('.model-region summary span:last-child').forEach(function (el) {
    var parts = el.textContent.split(' · ');
    el.textContent = parts[0].replace(' models', ' 3D models').replace(' model', ' 3D model');
  });
  document.querySelectorAll('.region-card dd').forEach(function (el) {
    var text = el.textContent;
    if (text.indexOf('summary row') !== -1) text = text.split(' · ')[0];
    text = text.replace('active PDBs', '3D models').replace('active PDB', '3D model');
    el.textContent = text;
  });
});
