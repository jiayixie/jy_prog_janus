set message = $argv[1]
git add -A
git commit -m "$message"
git push -u jy_prog_janus master # push the master barnch to jy_prog_janus repo.
