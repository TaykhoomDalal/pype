(() => {
  const root = document.documentElement;
  let savedTheme = null;
  try {
    savedTheme = localStorage.getItem("pype-theme");
  } catch {
    savedTheme = null;
  }
  if (savedTheme === "light" || savedTheme === "dark") {
    root.dataset.theme = savedTheme;
  }

  const currentTheme = () => {
    if (root.dataset.theme) {
      return root.dataset.theme;
    }
    return window.matchMedia("(prefers-color-scheme: dark)").matches
      ? "dark"
      : "light";
  };

  const updateButton = (button) => {
    const nextTheme = currentTheme() === "dark" ? "light" : "dark";
    button.textContent = nextTheme === "dark" ? "Dark mode" : "Light mode";
    button.setAttribute("aria-label", `Switch to ${nextTheme} mode`);
  };

  window.addEventListener("DOMContentLoaded", () => {
    document.querySelectorAll("[data-theme-toggle]").forEach((button) => {
      updateButton(button);
      button.addEventListener("click", () => {
        root.dataset.theme = currentTheme() === "dark" ? "light" : "dark";
        try {
          localStorage.setItem("pype-theme", root.dataset.theme);
        } catch {
          // The theme still applies for the current page.
        }
        document.querySelectorAll("[data-theme-toggle]").forEach(updateButton);
      });
    });
  });
})();
