(() => {
  const payloadNode = document.getElementById("ofp-report-payload");
  if (!payloadNode) {
    return;
  }

  const report = JSON.parse(payloadNode.textContent || "{}");
  const chainPanels = Array.isArray(report.chainPanels) ? report.chainPanels : [];
  if (!chainPanels.length) {
    return;
  }

  const reportPage = document.querySelector("[data-report-page]");
  const chainSelect = document.querySelector("[data-chain-select]");
  const activePanel = document.querySelector("[data-active-chain-panel]");
  const figureSheet = document.querySelector("[data-figure-sheet]");
  const summaryGrid = document.querySelector("[data-summary-grid]");
  const chainWarnings = document.querySelector("[data-chain-warnings]");
  const chainWarningsTitle = document.querySelector("[data-chain-warnings-title]");
  const chainWarningList = document.querySelector("[data-chain-warning-list]");
  const jobWarnings = document.querySelector("[data-job-warnings]");
  const jobWarningsTitle = document.querySelector("[data-job-warnings-title]");
  const jobWarningList = document.querySelector("[data-job-warning-list]");

  const panelByChain = new Map(chainPanels.map((panel) => [panel.referenceChain, panel]));

  const clearChildren = (node) => {
    if (!node) {
      return;
    }
    while (node.firstChild) {
      node.removeChild(node.firstChild);
    }
  };

  const buildSummaryItem = (item) => {
    const wrapper = document.createElement("div");
    wrapper.className = "ofp-summary-item summary-item";

    const labelRow = document.createElement("div");
    labelRow.className = "ofp-summary-label-row summary-label-row";

    const label = document.createElement("span");
    label.className = "ofp-summary-label summary-label";
    label.textContent = item.label;

    labelRow.append(label);

    if (item.tooltip) {
      const help = document.createElement("span");
      help.className = "ofp-summary-help summary-help";
      help.textContent = "?";
      help.title = item.tooltip;
      help.setAttribute("aria-label", item.tooltip);
      labelRow.append(help);
    }

    const value = document.createElement("span");
    value.className = "ofp-summary-value summary-value";
    value.textContent = item.value;

    wrapper.append(labelRow, value);
    return wrapper;
  };

  const renderSummary = (items) => {
    if (!summaryGrid) {
      return;
    }
    clearChildren(summaryGrid);
    if (!items.length) {
      return;
    }

    const alignedItems = items.slice(0, 8);
    const trailingItems = items.slice(8);

    const matrix = document.createElement("div");
    matrix.className = "ofp-summary-matrix summary-matrix";
    alignedItems.forEach((item) => {
      matrix.append(buildSummaryItem(item));
    });
    summaryGrid.append(matrix);

    if (trailingItems.length) {
      const trailing = document.createElement("div");
      trailing.className = "ofp-summary-trailing summary-trailing";
      trailingItems.forEach((item) => {
        trailing.append(buildSummaryItem(item));
      });
      summaryGrid.append(trailing);
    }
  };

  const renderWarnings = (card, titleNode, listNode, title, warnings) => {
    if (!card || !titleNode || !listNode) {
      return;
    }
    const items = Array.isArray(warnings) ? warnings : [];
    card.hidden = items.length === 0;
    if (items.length === 0) {
      clearChildren(listNode);
      return;
    }
    titleNode.textContent = items.length > 1 ? `${title} (${items.length})` : title;
    clearChildren(listNode);
    items.forEach((warning) => {
      const listItem = document.createElement("li");
      listItem.textContent = warning;
      listNode.append(listItem);
    });
    card.open = true;
  };

  const findTemplate = (chainId) => {
    const templates = document.querySelectorAll("template[data-chain-figure]");
    for (const template of templates) {
      if (template.dataset.chainFigure === chainId) {
        return template;
      }
    }
    return null;
  };

  const mountFigure = (chainId) => {
    if (!figureSheet) {
      return;
    }
    clearChildren(figureSheet);
    const template = findTemplate(chainId);
    if (!template) {
      return;
    }
    figureSheet.append(template.content.cloneNode(true));
  };

  const syncLocationHash = (chainId) => {
    const nextHash = `#chain-${chainId}`;
    try {
      window.history.replaceState(null, "", nextHash);
    } catch (_error) {
      window.location.hash = nextHash;
    }
  };

  const activateChain = (chainId, pushHash = true) => {
    const panel = panelByChain.get(chainId) || chainPanels[0];
    if (!panel) {
      return;
    }

    if (chainSelect) {
      chainSelect.value = panel.referenceChain;
    }
    if (activePanel) {
      activePanel.dataset.chainId = panel.referenceChain;
      activePanel.dataset.panelWidth = String(panel.panelWidth);
      activePanel.style.setProperty("--panel-width", `${panel.panelWidth}px`);
    }
    if (reportPage) {
      reportPage.style.setProperty("--active-panel-width", `${panel.panelWidth}px`);
    }

    renderSummary(panel.summaryItems || []);
    renderWarnings(chainWarnings, chainWarningsTitle, chainWarningList, "Current Chain Notes", panel.warnings || []);
    mountFigure(panel.referenceChain);

    if (pushHash) {
      syncLocationHash(panel.referenceChain);
    }
  };

  if (chainSelect) {
    clearChildren(chainSelect);
    chainPanels.forEach((panel) => {
      const option = document.createElement("option");
      option.value = panel.referenceChain;
      option.textContent = panel.chainLabel;
      chainSelect.append(option);
    });
    chainSelect.addEventListener("change", (event) => activateChain(event.target.value));
  }

  renderWarnings(jobWarnings, jobWarningsTitle, jobWarningList, "Job-level Notes", report.warnings || []);

  const requestedChain = window.location.hash.replace("#chain-", "");
  const initialChain = panelByChain.has(requestedChain) ? requestedChain : report.defaultReferenceChain;
  activateChain(initialChain || chainPanels[0].referenceChain, false);
})();
