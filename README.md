> <B>⚠️ This README is still under development</B>


<h1 align="center">
  <img src="https://github.com/CatarinaRRF/siRNA_Seeker/blob/main/static/img/siRNAseeker_preto.png" alt="logo">
</h1>

<p align="center">
    <img src="https://img.shields.io/github/last-commit/CatarinaRRF/siRNA_Seeker?color=informational&style=flat-square"
         alt="GitHub last commit">
    <a href="https://github.com/CatarinaRRF/siRNA_Seeker">
    <img src= http://img.shields.io/static/v1?label=STATUS&message=Beta&color=green&style=flat-square >
    </a>

</p>

<p align="center">
  <a href="#about">About</a> •
  <a href="#getting-started">Getting Started</a> •
  <a href="#software">Software</a> •
  <a href="#versioning">Versioning</a> •
  <a href="#license">License</a> •
  <a href="#credits">Credits</a>
</p>

### 🌌 Adjustments and Improvements

The project is still in development, and upcoming updates will focus on the following tasks:

- [ ] Documentation setup

# About 
<p align="justify"> The siRNA Seeker is a web app developed using the Django framework in Python to make available the siRNA design algorithm developed by the Nanogene research group to the scientific community. This web application offers a user-friendly and accessible interface for researchers to efficiently utilize the siRNA selection algorithm. Through this web app, users can interactively perform siRNA design analyses, customize parameters, and receive precise and relevant results. Django's use ensures robust and scalable implementation, meeting the demands of the scientific community and advancing biotechnology research. The goal is to make the web app accessible online, providing a centralized platform for knowledge sharing and collaboration among scientific community members interested in siRNA design. Thus, the following parts of the application have already been developed:

> ⚙️ Access more information about the algorithm developed [here](https://github.com/CatarinaRRF/IC_design_de_siRNA)

1. Landing Page - Main:
There is a banner with an image and a button to start siRNA design, which leads to the Form (where users input parameters).
Information about the algorithm, its importance for siRNA design, a contact option, and information about the research group are presented.
2. Form:
Users can input fasta file, copy and paste text, or put in the tag to search.
Users could choose the author, siRNA size, calculate Tm, ™ max, energy-free range, and minimum of compliance.
When clicking the submit button, users are directed to the Loading app.
3. Loading:
Show the expected time to complete the task and the progress.
When complete, users are directed to the Dashboard.
4. Dashboard:
Displays the gene being analyzed, the amount of selected sequences, and a comparative table for all possible sequences.
Text blocks indicate parameter failures, describing and classifying gravity.
5. Account:
Standard page for user login and registration. Accounts can save tree results at once.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To start this app in ypur local machine, you will need the following:

- [Django](https://docs.djangoproject.com/en/3.2/) - Web framework in Python
- [Celery](https://docs.celeryproject.org/en/stable/django/index.html) - For background task execution

If your operating system is Windows or MacOS, you will also need the following programs:
- [Redis](https://redis.io/documentation) - Messaging and cache system
- Linux emulator like [WSL](https://docs.microsoft.com/en-us/windows/wsl/install)

### Installation

A step-by-step series of examples that show how to get a development environment running

🛠️ <b>Clone this repository to your local environment</b>

<code>git clone https://github.com/CatarinaRRF/siRNA_Seeker.git</code>

🛠️ <b>Navigate to the project directory</b>

<code>cd siRNA_Seeker</code>

🛠️ <b>Install the necessary dependencies</b>

<code>pip install -r requirements.txt</code>

## Software

* [Django](https://www.djangoproject.com/) - The web framework used
* [Python](https://www.python.org/) - Used for algorithm production
* [Javascript]() - Used for chart production

## Versioning

The version control system used was [Git](https://git-scm.com/). 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Credits
* Federal University of Uberlândia (UFU) - Campus Patos de Minas
* Development team: Catarina RRF, Valdeir de Paula, and Matheus Souza
<img src="https://github.com/CatarinaRRF/Challenge-Alura-Cash-19-08-22/blob/974dd832c3980dd107a36a4b6906b616bb7b71f2/media/hr_line_redme.png" alt="logo">
<p align="center">

  <img src="https://github.com/CatarinaRRF/siRNA_Seeker/blob/main/static/img/logo_light.svg" alt="logo" width='150px'>
