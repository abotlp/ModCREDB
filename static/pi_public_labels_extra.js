(() => {
  const pairs = [
    ["Nearest Neighbor 50-70%", "Nearest Neighbor (70% - 40%)"],
    ["Nearest Neighbor 70-40%", "Nearest Neighbor (70% - 40%)"],
    ["Nearest Neighbor >70%", "Nearest Neighbor (>70%)"]
  ];
  const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  const nodes = [];
  while (walker.nextNode()) nodes.push(walker.currentNode);
  nodes.forEach(node => {
    let text = node.nodeValue;
    pairs.forEach(pair => { text = text.split(pair[0]).join(pair[1]); });
    node.nodeValue = text;
  });
})();
