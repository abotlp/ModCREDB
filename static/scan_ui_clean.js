document.addEventListener('DOMContentLoaded', function () {
  if (!location.pathname.startsWith('/scan')) return;

  var style = document.createElement('style');
  style.textContent = [
    '.scan-disclosure{background:#fff;border:1px solid var(--line);border-radius:10px;margin:18px 0;padding:0;overflow:hidden}',
    '.scan-disclosure>summary{cursor:pointer;display:flex;align-items:center;justify-content:space-between;gap:14px;list-style:none;padding:14px 16px;font-weight:750;background:#f7f9fc}',
    '.scan-disclosure>summary::-webkit-details-marker{display:none}',
    '.scan-disclosure>summary:before{content:"▸";display:inline-block;margin-right:8px;color:var(--accent);transition:transform .12s ease}',
    '.scan-disclosure[open]>summary:before{transform:rotate(90deg)}',
    '.scan-disclosure-body{padding:16px}',
    '.scan-disclosure .muted{font-weight:400}',
    '.scan-form .scan-disclosure{margin:10px 0}',
    '.scan-form .scan-disclosure>summary{background:#fbfcfe}',
    '.motif-picker:not([open]) .motif-picker-grid,.motif-picker:not([open]) form{display:none}'
  ].join('\n');
  document.head.appendChild(style);

  function makeSectionCollapsible(headingText, summaryText, noteText) {
    var headings = Array.from(document.querySelectorAll('h2,h3'));
    var heading = headings.find(function (h) { return h.textContent.trim().indexOf(headingText) !== -1; });
    if (!heading) return;
    var section = heading.closest('section') || heading.parentElement;
    if (!section || section.querySelector(':scope > details.scan-disclosure')) return;

    var details = document.createElement('details');
    details.className = 'scan-disclosure';
    var summary = document.createElement('summary');
    summary.innerHTML = '<span>' + summaryText + '</span>' + (noteText ? '<span class="muted">' + noteText + '</span>' : '');
    var body = document.createElement('div');
    body.className = 'scan-disclosure-body';

    var children = Array.from(section.childNodes);
    children.forEach(function (child) {
      if (child === heading) return;
      body.appendChild(child);
    });
    section.textContent = '';
    details.appendChild(summary);
    details.appendChild(body);
    section.appendChild(details);
  }

  makeSectionCollapsible('Filters and optional focused scan', 'Filters and optional focused scan', 'optional');

  document.querySelectorAll('details.motif-picker').forEach(function (details) {
    details.removeAttribute('open');
    details.classList.add('scan-disclosure');
    var summary = details.querySelector('summary');
    if (summary) {
      summary.innerHTML = '<span>Find individual motifs</span><span class="muted">optional</span>';
    }
    if (!details.querySelector('.scan-disclosure-body')) {
      var body = document.createElement('div');
      body.className = 'scan-disclosure-body';
      Array.from(details.childNodes).forEach(function (child) {
        if (child !== summary) body.appendChild(child);
      });
      details.appendChild(body);
    }
  });

  document.querySelectorAll('.motif-picker[open]').forEach(function (d) { d.removeAttribute('open'); });
});
