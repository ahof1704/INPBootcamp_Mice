<h1 align="center">VocalMat</h1>
<div align="center">
    <strong>INP Bootcamp 2020</strong>
</div>

<div align="center">
    <sub>This repository was built by Katie Ferguson, Clayton Barnes and Antonio Fonseca @ Neuroscience Department, Yale University.
</div>

<div align="center">
    <br />
    <!-- MATLAB version -->
    <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://img.shields.io/badge/MATLAB-2017a%7C2018b-blue.svg?style=flat-square"
      alt="MATLAB tested versions" />
    </a>
    <!-- LICENSE -->
    <a href="#">
    <img src="https://img.shields.io/badge/license-Apache%202.0-orange.svg?style=flat-square"
      alt="Apache License 2.0" />
    </a>
    <br />
</div>

## Table of Contents
- [Description](#description)
- [Getting Started](#getting-started)


## Description
> **This repository is organized according to the self-reported programming skill level of the INP students.**

<p align="justify"> The repository is divided into three main levels. The Beginner, Intermediary and Advanced.


<p align="justify"> <b>Beginner</b> [add description for Beginner programming level and what they will be working on]

<p align="justify"> <b>Intermediary</b> [add description for Intermediary programming level and what they will be working on]

<p align="justify"> <b>Advanced</b> [add description for Advanced programming level and what they will be working on]


## Getting Started

##### Cloning from the terminal (master)
```bash
$ git https://github.com/ahof1704/INPBootcamp_Mice.git
```

##### Using a Git client
You can use a Git client to clone our repository, we recommend GitHub's own client:
```
Download at: https://desktop.github.com
```

#### Directory Structure
- __Beginner:__ Directory with all files and scripts related to the VocalMat Identifier
- __Intermediary:__ Directory with all files and scripts related to the VocalMat Classifier
- __Advanced:__ Directory with audio files you want to analyze in the `audios` directory
- __Data:__ Directory with audio files you want to analyze in the `audios` directory
- __Plots:__ Directory with audio files you want to analyze in the `audios` directory

## Usage

#### `VocalMat` Output Files

<p align="justify">VocalMat outputs a directory with the same name as the audio file that was analyzed. Inside that directory there will be two directories (<i>All</i>, <i>All_axes</i>), and two Microsoft Excel (.xlsx) files. Inside <i>All_axes</i> you will find one image for each vocalization candidate detetcted with the resulting segmentation illusrated by sparsed blue dots. The raw original images are available inside <i>All</i>. The main Excel file has the same name of the audio file analyzed (<i>audio_file_name</i>.xlsx). This file contains information on each vocalization, such as start and end time, duration, frequency (minimum, mean and maximum), bandwidth, intensity (minimum, mean, maximum and corrected based on the backgroun), existence of harmonic components or distortions (noisy) and call type. The second excel file named as <i>audio_file_name</i>_DL.xlsx shows the probability distribution for each vocalization candidate across the different labels.

##### Software Requirements
- __MATLAB:__ versions 2017a through 2018b. For other versions refer to the [FAQ](#faq).
- __MATLAB Add-Ons:__
    - None

## License
<div>
    <a href="#">
    <img src="https://img.shields.io/badge/license-Apache%202.0-orange.svg?style=flat-square"
      alt="Apache License 2.0" />
    </a>
</div>

- **[Apache License 2.0](https://github.com/ahof1704/VocalMat/blob/VocalMat_RC/LICENSE)**
- Copyright 2019 © <a href="http://www.dietrich-lab.org" target="_blank">Dietrich Lab</a>.

<!-- version-control: 1.0 -->
