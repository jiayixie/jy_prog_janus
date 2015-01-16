# initiate a repo here, only use once when building new repo
set message = $argv[1]
git init
git config --global user.name "jiayixie"
git config --global user.email jiayi.seis@gmail.com
git add -A
git commit -m "$message" 
git remote add jy_prog_janus git@github.com:jiayixie/jy_prog_janus.git # use a name to represent the remote URL repo
git push -u jy_prog_janus master # push the master barnch to jy_prog_janus repo.

git pull jy_prog_janus master 
