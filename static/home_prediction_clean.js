document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('#examples tbody tr').forEach(function (row) {
    var cell = row.children[2];
    if (!cell) return;
    cell.querySelectorAll('.muted').forEach(function (el) { el.remove(); });
    cell.querySelectorAll('br').forEach(function (el) { el.remove(); });
  });
});
