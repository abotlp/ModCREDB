document.addEventListener('DOMContentLoaded', function () {
  if (location.pathname.indexOf('/scan') !== 0) return;
  document.querySelectorAll('details.motif-picker').forEach(function (box) {
    box.removeAttribute('open');
    var summary = box.querySelector('summary');
    if (summary) summary.setAttribute('title', 'Open optional motif search');
  });
});
