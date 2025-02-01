# MAE103

This repository is designed for MAE 103 Elementary Fluid Mechanics at UCLA.

[![Build Status](https://github.com/UCLAMAEThreads/MAE103/actions/workflows/CI.yml/badge.svg)](https://github.com/UCLAMAEThreads/MAE103/actions/workflows/CI.yml) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/UCLAMAEThreads/MAE103/HEAD)

You have two options for using the tools and notebooks in this repository. The first option is to run them in a cloud-based format, accessed by the Binder button above. In this format, you don't need to install anything locally, though you will need to download your notebooks to your own computer so that you can use them later. This process of downloading and uploading can be a bit annoying.

The other option, which gives you more freedom to work, is to run on your own computer. For this, you will need to follow these instructions to install jupyter and the julia programming language. To get started, please do the following:

1. Download Julia to your computer from https://julialang.org

2. Start a Julia session.

3. Type `]` to enter the Julia package management system. At this prompt, type `add IJulia`. This downloads `jupyter`, which is a browser-based computing environment. The 'ju' part stands for Julia; the 'py' part stands for Python. You can return to the Julia prompt by pressing backspace.

If you already have jupyter set up on your computer, start here:

1. Enter the Julia package management system using `]`, if you are not already in it (see above). Type `add https://github.com/UCLAMAEThreads/MAE103`. This will add the MAE103 package to your own computer. (In future sessions, you do not need to run this again. However, you should regularly run `update` to make sure you have the most up-to-date version of the package. I.e., type `] update` from the Julia prompt. Please make sure to start a new Julia session after you update, to force a new compilation of the package.)

2. Now press backspace to return to the Julia prompt. At this prompt, type `using MAE103`. This precompiles the package. Note: it will probably take several minutes.

3. After you have waited patiently for that to finish, type `MAE103.open_notebooks()`. This will put you into the main index of notebooks. You can open each one and run it.

4. If you wish to use any of these notebooks as a starting point for your own modifications, please use the "Make a Copy..." option from the main jupyter menu.


If you cannot get the notebooks to run on your own computer, you can run them on a remote server using the launch/binder button above. You will have to be careful to download your saved notebooks to your own computer, since you generally cannot save your work.

