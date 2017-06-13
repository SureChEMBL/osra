"master" branch is for upstream
"deb" branch for debianization and SureChEMBL specific changes

If you want to setup svn-remote for fetching the latest changes from
the original SVN repository, add the following lines in .git/config:

[svn-remote "svn"]
        url = svn://svn.code.sf.net/p/osra/code
        fetch = trunk:refs/remotes/svn/trunk
        branches = branches/*:refs/svn/upsteram/*
        tags = tags/*:refs/remotes/svn/tags/*
[svn]
        authorsfile = ./authors.txt

After that "git branch -r" will show "svn/<smth>" remote branches.

For building debian packages do the following:
$ rm -rf debian
$ ./configure
$ cp -r package/linux/debian ./
$ debuild -b
