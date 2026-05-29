# Burhan Beycan - Personal CV Website

A professional, responsive CV website showcasing the research, AI, and academic achievements of Burhan Beycan, a multidisciplinary researcher with experience across applied AI, computer vision, polymer chemistry, and materials science.

## 🎯 Features

- **Responsive Design**: Optimized for desktop, tablet, and mobile devices
- **Modern UI/UX**: Clean, professional design with smooth animations
- **Interactive Portfolio**: Filterable research project showcase
- **GitHub Project Integration**: Dedicated project map for AI, computer vision, scientific AI, and research tools
- **Publications Section**: Comprehensive list of academic publications
- **Contact Form**: Functional contact form with validation
- **Performance Optimized**: Fast loading with optimized assets

## 🧭 Career Positioning

Recommended positioning for the website and GitHub portfolio:

> Applied AI / Machine Learning Researcher with PhD-level scientific R&D expertise

The GitHub integration is designed to avoid limiting the profile to only AI-for-chemistry or AI-for-materials. It highlights three career signals:

1. **General applied AI**: RAG, evaluation, deployment, APIs, and product-ready ML applications.
2. **Computer vision**: detection, tracking, image analytics, and real-time model workflows.
3. **Scientific AI**: materials informatics, experimental design, Bayesian optimization, and AI for R&D.

## 🔗 GitHub Integration

The dedicated GitHub project page is available at:

```text
github-projects.html
```

It loads project data from:

```text
assets/data/github-projects.json
```

The page renders cards dynamically using:

```text
assets/js/github-projects.js
assets/css/github-projects.css
```

The project data is intentionally separated from the HTML so the portfolio can be updated by editing one JSON file instead of manually editing project cards.

## 🛠️ Technologies Used

- **HTML5**: Semantic markup and accessibility
- **CSS3**: Modern styling with animations and responsive design
- **JavaScript**: Interactive functionality, GitHub project rendering, and form handling
- **Bootstrap 5**: Responsive grid system and components
- **GitHub REST API**: Repository metadata enrichment on the GitHub projects page
- **Python**: Project data validation script
- **GitHub Actions**: Automated validation for the GitHub project data and integration assets
- **AOS**: Animate On Scroll library
- **Typed.js**: Dynamic typing animations
- **Isotope**: Portfolio filtering
- **GLightbox**: Image lightbox functionality
- **Swiper**: Touch slider for portfolio details

## 📁 Project Structure

```text
cv-website/
├── index.html                              # Main homepage
├── github-projects.html                    # GitHub-integrated AI project map
├── portfolio-details.html                  # Portfolio detail page
├── assets/
│   ├── css/
│   │   ├── style.css                       # Main stylesheet
│   │   └── github-projects.css             # GitHub projects page styles
│   ├── js/
│   │   ├── main.js                         # Main JavaScript file
│   │   └── github-projects.js              # Dynamic GitHub project renderer
│   ├── data/
│   │   └── github-projects.json            # Project portfolio data source
│   ├── img/
│   │   ├── profile-img.jpg                 # Profile photo
│   │   ├── hero-bg.jpg                     # Hero background
│   │   ├── favicon.png                     # Website favicon
│   │   └── portfolio/                      # Portfolio images
│   └── cv/
│       └── CV_BBEYCAN_CENG.pdf             # Downloadable CV
├── scripts/
│   └── validate-github-projects.py         # JSON validator for project data
├── .github/
│   └── workflows/
│       └── validate-github-projects.yml    # GitHub Actions validation workflow
├── docs/
│   ├── GITHUB_INTEGRATION.md               # Maintenance and editing guide
│   └── PROJECT_PORTFOLIO_STRATEGY.md       # Career and project roadmap
├── forms/
│   └── contact.php                         # Contact form handler
└── README.md                               # This file
```

## ✅ Maintaining the GitHub Project Data

After editing `assets/data/github-projects.json`, run:

```bash
python scripts/validate-github-projects.py assets/data/github-projects.json
```

The validator checks required fields, duplicate IDs, valid categories, valid statuses, URL format, and list fields. The same check runs automatically in GitHub Actions for relevant pull requests and pushes to `main`.

## 📱 Browser Support

- Chrome (latest)
- Firefox (latest)
- Safari (latest)
- Edge (latest)
- Mobile browsers (iOS Safari, Chrome Mobile)

## 🎨 Design Credits

- Template inspired by [BootstrapMade](https://bootstrapmade.com/)
- Icons from [Bootstrap Icons](https://icons.getbootstrap.com/) and [Boxicons](https://boxicons.com/)
- Fonts from [Google Fonts](https://fonts.google.com/)

## 📄 License

This project is open source and available under the [MIT License](LICENSE).

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](../../issues).

## 📞 Support

If you have any questions or need help with customization, please feel free to reach out:

- Email: burhanbeycan@hotmail.com
- LinkedIn: [linkedin.com/in/burhanbeycan](https://www.linkedin.com/in/burhanbeycan)

## 🌟 Acknowledgments

- Thanks to the open-source community for the amazing libraries and tools
- Special thanks to BootstrapMade for the design inspiration
- Research institutions and collaborators for their support

---

**Made with ❤️ by Burhan Beycan**
