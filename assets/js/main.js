/**
* Template Name: iPortfolio - v3.10.0
* Template URL: https://bootstrapmade.com/iportfolio-bootstrap-portfolio-websites-template/
* Author: BootstrapMade.com
* License: https://bootstrapmade.com/license/
* Modified for Burhan Beycan CV Website
*/

(function() {
  "use strict";

  /**
   * Easy selector helper function
   */
  const select = (el, all = false) => {
    el = el.trim()
    if (all) {
      return [...document.querySelectorAll(el)]
    } else {
      return document.querySelector(el)
    }
  }

  /**
   * Easy event listener function
   */
  const on = (type, el, listener, all = false) => {
    let selectEl = select(el, all)
    if (selectEl) {
      if (all) {
        selectEl.forEach(e => e.addEventListener(type, listener))
      } else {
        selectEl.addEventListener(type, listener)
      }
    }
  }

  /**
   * Easy on scroll event listener 
   */
  const onscroll = (el, listener) => {
    el.addEventListener('scroll', listener)
  }

  /**
   * Navbar links active state on scroll
   */
  let navbarlinks = select('#navbar .scrollto', true)
  const navbarlinksActive = () => {
    let position = window.scrollY + 200
    navbarlinks.forEach(navbarlink => {
      if (!navbarlink.hash) return
      let section = select(navbarlink.hash)
      if (!section) return
      if (position >= section.offsetTop && position <= (section.offsetTop + section.offsetHeight)) {
        navbarlink.classList.add('active')
      } else {
        navbarlink.classList.remove('active')
      }
    })
  }
  window.addEventListener('load', navbarlinksActive)
  onscroll(document, navbarlinksActive)

  /**
   * Scrolls to an element with header offset
   */
  const scrollto = (el) => {
    let elementPos = select(el).offsetTop
    window.scrollTo({
      top: elementPos,
      behavior: 'smooth'
    })
  }

  /**
   * Back to top button
   */
  let backtotop = select('.back-to-top')
  if (backtotop) {
    const toggleBacktotop = () => {
      if (window.scrollY > 100) {
        backtotop.classList.add('active')
      } else {
        backtotop.classList.remove('active')
      }
    }
    window.addEventListener('load', toggleBacktotop)
    onscroll(document, toggleBacktotop)
  }

  /**
   * Mobile nav toggle
   */
  on('click', '.mobile-nav-toggle', function(e) {
    select('body').classList.toggle('mobile-nav-active')
    this.classList.toggle('bi-list')
    this.classList.toggle('bi-x')
  })

  /**
   * Scrool with ofset on links with a class name .scrollto
   */
  on('click', '.scrollto', function(e) {
    if (select(this.hash)) {
      e.preventDefault()

      let body = select('body')
      if (body.classList.contains('mobile-nav-active')) {
        body.classList.remove('mobile-nav-active')
        let navbarToggle = select('.mobile-nav-toggle')
        navbarToggle.classList.toggle('bi-list')
        navbarToggle.classList.toggle('bi-x')
      }
      scrollto(this.hash)
    }
  }, true)

  /**
   * Scroll with ofset on page load with hash links in the url
   */
  window.addEventListener('load', () => {
    if (window.location.hash) {
      if (select(window.location.hash)) {
        scrollto(window.location.hash)
      }
    }
  });

  /**
   * Hero typed text effect
   */
  const typed = select('.typed')
  if (typed) {
    let typed_strings = typed.getAttribute('data-typed-items')
    typed_strings = typed_strings.split(',')
    new Typed('.typed', {
      strings: typed_strings,
      loop: true,
      typeSpeed: 100,
      backSpeed: 50,
      backDelay: 2000
    });
  }

  /**
   * Skills animation
   */
  let skillsSection = select('#skills');
  if (skillsSection) {
    let skillsIsAnimated = false;
    
    const animateSkills = () => {
      if (!skillsIsAnimated) {
        let skillsContent = select('.skills-content');
        if (skillsContent) {
          let progressBars = select('.progress .progress-bar', true);
          progressBars.forEach((progressBar) => {
            let progress = progressBar.getAttribute('aria-valuenow');
            progressBar.style.width = progress + '%';
          });
          skillsIsAnimated = true;
        }
      }
    };

    // Check if skills section is in viewport
    const checkSkillsInView = () => {
      let skillsTop = skillsSection.offsetTop;
      let skillsHeight = skillsSection.offsetHeight;
      let windowTop = window.pageYOffset;
      let windowHeight = window.innerHeight;

      if (windowTop + windowHeight >= skillsTop && windowTop <= skillsTop + skillsHeight) {
        animateSkills();
      }
    };

    window.addEventListener('scroll', checkSkillsInView);
    window.addEventListener('load', checkSkillsInView);
  }

  /**
   * Porfolio isotope and filter
   */
  window.addEventListener('load', () => {
    let portfolioContainer = select('.portfolio-container');
    if (portfolioContainer) {
      let portfolioIsotope = new Isotope(portfolioContainer, {
        itemSelector: '.portfolio-item'
      });

      let portfolioFilters = select('#portfolio-flters li', true);

      on('click', '#portfolio-flters li', function(e) {
        e.preventDefault();
        portfolioFilters.forEach(function(el) {
          el.classList.remove('filter-active');
        });
        this.classList.add('filter-active');

        portfolioIsotope.arrange({
          filter: this.getAttribute('data-filter')
        });
        portfolioIsotope.on('arrangeComplete', function() {
          AOS.refresh()
        });
      }, true);
    }

  });

  /**
   * Initiate portfolio lightbox 
   */
  const portfolioLightbox = GLightbox({
    selector: '.portfolio-lightbox'
  });

  /**
   * Portfolio details slider
   */
  new Swiper('.portfolio-details-slider', {
    speed: 400,
    loop: true,
    autoplay: {
      delay: 5000,
      disableOnInteraction: false
    },
    pagination: {
      el: '.swiper-pagination',
      type: 'bullets',
      clickable: true
    }
  });

  /**
   * Animation on scroll
   */
  window.addEventListener('load', () => {
    AOS.init({
      duration: 1000,
      easing: 'ease-in-out',
      once: true,
      mirror: false
    })
  });

  /**
   * Initiate Pure Counter 
   */
  new PureCounter();

  /**
   * Contact form validation and submission
   */
  const contactForm = select('.php-email-form');
  if (contactForm) {
    contactForm.addEventListener('submit', function(e) {
      e.preventDefault();
      
      let thisForm = this;
      let action = thisForm.getAttribute('action');
      let recaptcha = thisForm.getAttribute('data-recaptcha-site-key');
      
      if (!action) {
        displayError(thisForm, 'The form action property is not set!');
        return;
      }
      
      thisForm.querySelector('.loading').classList.add('d-block');
      thisForm.querySelector('.error-message').classList.remove('d-block');
      thisForm.querySelector('.sent-message').classList.remove('d-block');
      
      let formData = new FormData(thisForm);
      
      if (recaptcha) {
        if (typeof grecaptcha !== "undefined") {
          grecaptcha.ready(function() {
            try {
              grecaptcha.execute(recaptcha, {action: 'php_email_form_submit'})
              .then(token => {
                formData.set('recaptcha-response', token);
                php_email_form_submit(thisForm, action, formData);
              })
            } catch(error) {
              displayError(thisForm, error);
            }
          });
        } else {
          displayError(thisForm, 'The reCaptcha javascript API url is not loaded!')
        }
      } else {
        php_email_form_submit(thisForm, action, formData);
      }
    });
  }

  function php_email_form_submit(thisForm, action, formData) {
    fetch(action, {
      method: 'POST',
      body: formData,
      headers: {'X-Requested-With': 'XMLHttpRequest'}
    })
    .then(response => {
      if (response.ok) {
        return response.text();
      } else {
        throw new Error(`${response.status} ${response.statusText} ${response.url}`);
      }
    })
    .then(data => {
      thisForm.querySelector('.loading').classList.remove('d-block');
      if (data.trim() == 'OK') {
        thisForm.querySelector('.sent-message').classList.add('d-block');
        thisForm.reset();
      } else {
        throw new Error(data ? data : 'Form submission failed and no error message returned from: ' + action);
      }
    })
    .catch((error) => {
      displayError(thisForm, error);
    });
  }

  function displayError(thisForm, error) {
    thisForm.querySelector('.loading').classList.remove('d-block');
    thisForm.querySelector('.error-message').innerHTML = error;
    thisForm.querySelector('.error-message').classList.add('d-block');
  }

  /**
   * Custom enhancements
   */

  // Smooth reveal animations for sections
  const observerOptions = {
    threshold: 0.1,
    rootMargin: '0px 0px -50px 0px'
  };

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        entry.target.classList.add('loaded');
      }
    });
  }, observerOptions);

  // Observe all sections for animation
  window.addEventListener('load', () => {
    const sections = select('section', true);
    sections.forEach(section => {
      section.classList.add('loading-animation');
      observer.observe(section);
    });
  });

  // Enhanced typing effect for hero section
  const heroTyped = select('#hero .typed');
  if (heroTyped) {
    // Add cursor blinking effect
    const style = document.createElement('style');
    style.textContent = `
      .typed-cursor {
        opacity: 1;
        animation: typedjsBlink 0.7s infinite;
      }
      @keyframes typedjsBlink {
        50% { opacity: 0; }
      }
    `;
    document.head.appendChild(style);
  }

  // Dynamic year update in footer
  const currentYear = new Date().getFullYear();
  const copyrightElement = select('#footer .copyright');
  if (copyrightElement && currentYear > 2025) {
    copyrightElement.innerHTML = copyrightElement.innerHTML.replace('2025', currentYear);
  }

  // Enhanced portfolio filtering with animation
  const portfolioItems = select('.portfolio-item', true);
  if (portfolioItems.length > 0) {
    portfolioItems.forEach((item, index) => {
      item.style.animationDelay = `${index * 0.1}s`;
    });
  }

  // Preload critical images
  const criticalImages = [
    'assets/img/profile-img.jpg',
    'assets/img/hero-bg.jpg'
  ];

  criticalImages.forEach(src => {
    const img = new Image();
    img.src = src;
  });

  // Add loading states for better UX
  window.addEventListener('load', () => {
    document.body.classList.add('loaded');
    
    // Remove loading overlay if exists
    const loadingOverlay = select('.loading-overlay');
    if (loadingOverlay) {
      setTimeout(() => {
        loadingOverlay.style.opacity = '0';
        setTimeout(() => {
          loadingOverlay.remove();
        }, 300);
      }, 500);
    }
  });

  // Enhanced scroll behavior for mobile
  let ticking = false;
  
  function updateScrollPosition() {
    navbarlinksActive();
    
    // Update progress indicator if exists
    const progressBar = select('.scroll-progress');
    if (progressBar) {
      const scrolled = (window.scrollY / (document.documentElement.scrollHeight - window.innerHeight)) * 100;
      progressBar.style.width = scrolled + '%';
    }
    
    ticking = false;
  }

  function requestTick() {
    if (!ticking) {
      requestAnimationFrame(updateScrollPosition);
      ticking = true;
    }
  }

  window.addEventListener('scroll', requestTick);

  // Add keyboard navigation support
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
      // Close mobile menu if open
      const body = select('body');
      if (body.classList.contains('mobile-nav-active')) {
        body.classList.remove('mobile-nav-active');
        const navbarToggle = select('.mobile-nav-toggle');
        navbarToggle.classList.add('bi-list');
        navbarToggle.classList.remove('bi-x');
      }
    }
  });

  // Enhanced form validation
  const formInputs = select('.php-email-form input, .php-email-form textarea', true);
  if (formInputs.length > 0) {
    formInputs.forEach(input => {
      input.addEventListener('blur', function() {
        validateField(this);
      });
      
      input.addEventListener('input', function() {
        if (this.classList.contains('is-invalid')) {
          validateField(this);
        }
      });
    });
  }

  function validateField(field) {
    const value = field.value.trim();
    const type = field.type;
    let isValid = true;
    let errorMessage = '';

    if (field.hasAttribute('required') && !value) {
      isValid = false;
      errorMessage = 'This field is required.';
    } else if (type === 'email' && value) {
      const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
      if (!emailRegex.test(value)) {
        isValid = false;
        errorMessage = 'Please enter a valid email address.';
      }
    }

    // Update field styling
    if (isValid) {
      field.classList.remove('is-invalid');
      field.classList.add('is-valid');
    } else {
      field.classList.remove('is-valid');
      field.classList.add('is-invalid');
    }

    // Show/hide error message
    let errorElement = field.parentNode.querySelector('.field-error');
    if (!errorElement && !isValid) {
      errorElement = document.createElement('div');
      errorElement.className = 'field-error text-danger small mt-1';
      field.parentNode.appendChild(errorElement);
    }
    
    if (errorElement) {
      errorElement.textContent = errorMessage;
      errorElement.style.display = isValid ? 'none' : 'block';
    }

    return isValid;
  }

})();
