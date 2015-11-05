# DirectionalDataSlides

First week's class notes.

The source file in all cases is the emacs
[org-mode](http://orgmode.org/) file.  Org can export directly to html
and markdown, and (via latex) to pdf (including beamer presentations).
These files are included in the repository.

Interestingly, it's not easy to view the html files on github.  One
work around is to use htmlpreview.github.com.  Just add the url of the
html file to htmlpreview.github.com/?, e.g.,

[http://htmlpreview.github.io/?https://github.com/presnell/DirectionalDataSlides/blob/master/directionalData.html]

Another interesting tool is [pandoc](http://pandoc.org/), which can
convert between many document formats.  In particular, pandoc will
convert emacs org-mode files to html, markdown, and pdf (again via
latex).  So for example, one can run the command
```sh
pandoc directionalData.org -o directionalData-pandoc.html
```
to produce html output, and
```sh
pandoc directionalData.org -o directionalData-pandoc.pdf
```
to produce pdf output.  Of course the results are not as good as can
be achieved using org-mode's native exporting facilities, but that
does not diminish pandoc's usefulness in a wide variety of contexts.

