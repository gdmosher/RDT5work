git checkout gh-pages
rsync -avh --stats --progress --delete-after --exclude .git/ ~/Dropbox/Websites/manuals/_site/ ~/Dropbox/Websites/manuals/
git add -A :/
git commit -am "some edits - see master commits for details"; git push -u origin gh-pages
git checkout master

