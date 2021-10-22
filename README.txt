+***Per runnare i jupiter e generare le figure è prima necessario generare i file contenenti le simulazioni necessarie tramite mol2.py***

I nome dei file del tipo fig11 si riferiscono alle figure già originariamente presenti nella tesi di Cavallini. Il notebook Fig_11 contiene dunque uno studio del caso riportato in Fig11 da Cavallini. Oltre alle figure presenti nella tesi, ne sono generate altre che non sono state scelte. Tutte queste figure sono nella cartella immagini, mentre quelle usate nella tesi sono nella cartella col codice latex.

Si consiglia sempre di eseguire "run all cells" se qualcosa sembra strano, in quanto le celle dei notebook non sono autoconsistenti. Anche una volta eseguito l'intero notebook, a volte, eseguire nuovamenteb una singola cella può dare un risultato scorretto se non è ripetuta l'esecuzione delle precedenti. 

L'unica modifica al programma mol2.py è l'inserimento del parametro g2, che consente di assegnare alla seconda particella un coefficiente di smorzamento diverso dalla prima. Se non specificato, entrambe le particelle avranno lo stesso gamma, g nel codice. 

La cartella "filmino" consente di generare la gif della molecola usata nella presentazione.



