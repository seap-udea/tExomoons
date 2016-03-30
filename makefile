# Para agregar archivos: git add <archivo>

push:
	-git commit -am "New changes"
	-git push origin master

pull:
	-git pull

clean:
	-find . -name "*~" -exec rm {} \;

