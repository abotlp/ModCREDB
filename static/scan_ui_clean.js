document.addEventListener('DOMContentLoaded', function () {
  if (location.pathname.indexOf('/scan') !== 0) return;

  var focus = document.querySelector('details.scan-focus');
  if (focus) {
    var summary = focus.querySelector('summary');
    if (summary) {
      summary.innerHTML = '<span>Personalized TF scan</span><span class="muted">optional</span>';
    }
    var label = Array.from(focus.querySelectorAll('label')).find(function (el) {
      return el.textContent.trim().indexOf('Selected motifs') === 0;
    });
    if (label) {
      var note = label.nextElementSibling;
      var motifBox = document.createElement('details');
      motifBox.className = 'scan-focus';
      if (label.querySelector('textarea') && label.querySelector('textarea').value.trim()) motifBox.setAttribute('open', '');
      motifBox.innerHTML = '<summary><span>Specific motif scan</span><span class="muted">optional</span></summary><div class="scan-focus-body"></div>';
      var body = motifBox.querySelector('.scan-focus-body');
      body.appendChild(label);
      if (note && note.classList && note.classList.contains('muted')) {
        note.textContent = 'Optional. Leave blank for the global generated-PWM scan.';
        body.appendChild(note);
      }
      focus.parentNode.insertBefore(motifBox, focus.nextSibling);
    }
  }

  document.querySelectorAll('details.motif-picker').forEach(function (box) {
    box.removeAttribute('open');
    var summary = box.querySelector('summary');
    if (summary) {
      summary.innerHTML = '<span>Find individual motifs</span><span class="muted">optional</span>';
      summary.setAttribute('title', 'Open optional motif search');
    }
  });
});
