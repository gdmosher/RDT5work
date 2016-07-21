---
title: Make the UCR Documentation Theme Your Theme!
sidebar: mydoc_sidebar
permalink: /README/
type: homepage
---
## UCR Documentation Theme

The UCR Documentation Theme is making reproducible research more accessible for all labs on campus. Simple to setup and maintain, easy to publish at github-pages, this theme supports advanced data publishing technology including R, R-Markdown and knitr. Development of the theme is supervised by [Thomas Girke](mailto:thomas.girke@ucr.edu) at UCR's IGGB/biocluster? The theme is supported by [Gordon David Mosher](mailto:gmosh001@ucr.edu) and for the time being, includes free onsite consultation and training to help you adopt the theme, import your data, and get started with markdown.

This theme is based on Tom Johnson's Jekyll Documentation theme version 5.0

## Jekyll Documentation Theme

The base theme continues to evolve here: [http://idratherbewriting.com/documentation-theme-jekyll/](http://idratherbewriting.com/documentation-theme-jekyll/) The instructions for the base theme 5.0 can be found by invoking "jekyll serve" in the root of the ORIG branch of this repository.



## Getting started

These instructions are specific to linux, but can be adapted to Windows or OSX.

Adding your content is one of the last steps in this tutorial. Follow along in order to get a top down view of the complete process.

### Prerequisites

To use this tool chain you will need to have already installed [git](https://help.github.com/articles/set-up-git/) and [jekyll](https://jekyllrb.com/docs/installation/). \\
If you are lucky, and using linux, the following commands will take care of that for you:

```
sudo -i
apt-get install git
git config --global user.name "YOUR NAME"
git config --global user.email "YOUR EMAIL ADDRESS"

gem install jekyll  #or apt-get install jekyll
```
To use the md2jekyllDT5 script to fold your files into sections you will also need [R](https://www.r-project.org/) and optionally [RStudio](https://www.rstudio.com/home/).\\
First, get R

```
sudo apt-get install r-base #installs R
```

Now [install RStudio using an installer](https://www.rstudio.com/products/rstudio/download/)
RStudio may ask for libraries/packages. Install them as needed at that time.
Additional help from [CS UCLA](http://web.cs.ucla.edu/~gulzar/rstudio/)

[Dropbox](https://www.dropbox.com/install) is recommended for multi-computer access and backup features.

### Clone or download the template

To facilitate use this theme as a template we have created a striped down configuration in the "gh-pages" branch. \\
To get your copy

- Note: Keeping your repo in a Dropbox folder protects your files and makes them accessible from your other computers!

- Open your shell and navigate to a folder on your computer where you want to create it,  
    then enter these commands to clone it and checkout the branch called "gh-pages":  
  
```
## Using this method the .git folder will exist although it may be hidden

git clone https://github.com/gdmosher/RDT5.git RDT5
cd RDT5
git checkout gh-pages #to open the gh-pages branch as a template
## -OR - git checkout ORIG-JDT5 #if you want to deviate from the tutorial
##                               to study Tom Johnson version 5.0 release
```

- All shell commands are assumed to be executed in the root of the repo unless otherwise noted.
- The root folder can be renamed at anytime, but I recommend leaving it until you are comfortable editing the _config.yml file.
    Important the baseurl: in _config.yml (described below), must match the repo name.

>If you want to start your repo fresh (Note: Don't do this unless you understand git and really want to get rid of all the history and other branches! You will lose all but the current branch, and it will become master.)
>
>-   delete the .git folder, then  
-   git init  
>
>To record this fresh state, you will want to
>
>-   git add .  
-   git commit -m "Initial commit of clean RDT5"
-   END of git init steps  


As you work in git be sure to stash or commit your changes before changing branches so your work doesn't get mixed up.
If you are new to git, this is a great time to create a branch for your work:

```
##   You can skip this step too, because the template is already in the gh-pages branch
##git checkout -b gh-pages  
##   Note: The name of the branch doesn't matter for now,
##   but some of the publishing scripts discussed below expect the name gh-pages
##   which is required by github to host your pages.
```
    
The site can be generated now, local and empty, if you want to test it.

```
jekyll serve
##  jekyll will generate an url for you to paste into your browser's address bar
```

Before you can push it to github you will need to create an empty repo at github.

-   at github.com login as yourself or create a free account
-   follow their instructions to create an empty repo
    -   do not add a README, License or .gitignore at this time
-   done

To connect your clone of RDT5 to your new repo at github

```
git remote add origin https://github.com/YourGithubID/YourRepoName.git
git remote -v    #to verify
```
-   when you are ready to push,  study the available publishing scripts [buildAll.sh](https://github.com/gdmosher/RDT5/blob/gh-pages/buildAll.sh) and [pushSite.sh](https://github.com/gdmosher/RDT5/blob/gh-pages/pushSite.sh)

-   then use one of the scripts or the command:

```
git push -u origin   #to push the branch you are on
```

### Architecture
Brief introduction to the architecture of the repository. The SIDEBAR is one of the hottest features of this theme. Multiple sidebars can be supported by one repo/site. The term PRODUCT refers to the content of one sidebar. Each product should be stored in it's own folder in the root of the site (e.g. /RDT5/mydoc). The default product name "mydoc" already has an empty folder ready for your content. If you change the name of the mydoc folder or want to add another product folder, you will need to learn to connect them by configuring several of the .yml and .html files in their various folders. In addition to populating the mydoc folder you will want to provide an index.md or index.html file in the root which will serve as a landing page. (Note: The folders "Product 1" and "Product 2" are already connected too.) (Note: jekyll won't try to generate pages in  folders with names that start with underscores.)
    
### Configuration
The template was designed so that no configuration is required to generate the empty site. But there are many placeholders for you to enter your project specific names and titles. The most important configuration files are listed in this section and the following shell command can help you locate them as well:  

```
grep -ri "your" *  
```

You can use your favorite command line or graphical tools to carefully edit the .yml and .html files to substitute your data. You shouldn't need it, but if you want more instructions on how to use the configuration files to change the behavior of the site, refer to Tom Johnson's version 5.0 documentation as described at the top of this page. Here is a list of files that you will want to understand first:  
    
-   [_config.yml](https://github.com/gdmosher/RDT5/blob/gh-pages/_config.yml)   
-   [_data/topnav.yml](https://github.com/gdmosher/RDT5/blob/gh-pages/_data/topnav.yml)  
-   [_data/sidebars/mydoc_sidebar.yml](https://github.com/gdmosher/RDT5/blob/gh-pages/_data/sidebars/mydoc_sidebar.yml)  
Note: For reference, these links point to the template online. To make changes you will edit in your copy of the repository.

if you want to create a new product folder, beyond the three already provided or with a name you like better, you will need to:

-   create the new folder in the root of the repo to hold your new markdown pages,
-   create a new prodName_sidebar.yml in _data/sidebars,
-   hook it up in [_includes/custom/sidebarconfigs.html](https://github.com/gdmosher/RDT5/blob/gh-pages/_includes/custom/sidebarconfigs.html),
-   and possibly configure some defaults for it in _config.yml

### Automation
Several shell scripts are provided in the root to help you learn to use the system and make some repetitive steps easier. It's always a good idea to make a backup before attempting new automation steps.

The command "jekyll build" will generate the site into the folder called _site. "jekyll serve" will build the site and serve it on the localhost [http://127.0.0.1:4005/RDT5/](http://127.0.0.1:4005/RDT5/) where RDT5 is the name of the repo folder and is the baseurl: in the _config.yml configuration file.  

Next we have three scripts that will push your site out to the web. These are buildAll.sh, pushSite.sh, and publish.sh which are discussed in the section [Publish It](#publish-it). Another level of automation called md2jekyllDT5 that folds your documents into sections will be discussed later in the section [Adding Your Content](#adding-your-content).

### Frontmatter
Documentation Theme frontmatter must be added to each .md file so that the theme can generate the links from the sidebar to the pages. If you use the md2jekyllDT5 script described in the next section, it will generate the sidebar and permalink frontmatter for you. The title is mandatory and it is your responsibility to place it in your .md files. Keywords and comments are optional, but you should use them as place holders. Currently the md2jekyllDT5 script inserts the sidebar/permalink frontmatter at lines 5-6  in your .md file and does not check that it falls inside the "- - -" fences. The frontmatter needs to be placed at the very beginning of the file and looks like this:  
    \-\-\-  
    title: Your Title  
    keywords:  
    comments:  
    sidebar: mydoc_sidebar  
    permalink: /folder/page/  
    \-\-\-  
    
If you have some files that don't use the md2jekyllDT5 script, you should add the sidebar and permalink frontmatter manually to your document and also to the mydoc_sidebar.yml. Permalinks are case sensitive and should generally contain the foldername and pagename with slashes before and after as shown above. 

### Adding Your Content
Jekyll expects your content in markdown (.md) format and the Documentation Theme expects it in the folder called mydocs, so just drop your content in there, add the frontmatter to the pages and the mydoc_sidebar.yml, and then type jekyll serve to generate the site. Done!

But there is more automation available if you want it. We copy in a folder called \_vignettes that contains many filename.knit.md files that are created by kniting in R with the option keep intermediate (or is it don't cleanup?). These files are submitted to the [md2jekyllDT5.R](https://github.com/gdmosher/RDT5/blob/gh-pages/_vignettes/md2jekyllDT5.R) script, also located in the _vignettes folder, to be folded into sections the output of which goes into the mydoc folder with their corresponding images. These *.knit.md files are each kept in their own folder. When starting the md2jekyllDT5 script you must pass in the several parameters. Find the parameters and an example of their use in the last few lines of the md2jekyllDT5.R script itself. I run the script in RStudio, but it can just as easily be run from the R command line.

### Publish It
Finally, we have three scripts that will push your site out to the web. These are: [buildAll.sh](https://github.com/gdmosher/RDT5/blob/gh-pages/buildAll.sh),   [pushSite.sh](https://github.com/gdmosher/RDT5/blob/gh-pages/pushSite.sh),  and  [publish.sh](https://github.com/gdmosher/RDT5/blob/gh-pages/publish.sh)

Please open the files and understand the commands inside before running them, because they may make significant changes to your repository. Determine which format and destination is best for you and then stick with it because they are not specifically designed to be mixed.

