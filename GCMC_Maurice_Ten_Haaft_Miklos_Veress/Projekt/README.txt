Für Windows:
--------------------------------------------
2.        Quellcode
--------------------------------------------
-In gcmc_vs.cpp sind alle Methoden mit der main methode enthalten
-gcmc_plotting enthält den python-skript zu der Generierung unserer Plots

------------------------------------------------
1. Starten einer Simulation ohne Kompilierung
------------------------------------------------
-Simulationen können... 
	...sequenziell über run_simulations.bat (Batch file) gestartet werden
	...parallel über run_simulations_parallel.bat gestartet werden
-die Batch files führen GCMC_VS.exe aus
-Es sind jeweils Beispiele eingetragen, die Parameter haben im batch-file die folgende Reihenfolge:
<.exe> <output_path + output_filename> <z> <thermalisation_steps> <interval> <total_measurements> <seed>

-Die Gesamtzahl an Iterationen ergibt sich dann aus thermalistation_steps + interval * total_measurements
-Seed ist für den Zufallsgenerator
-die Ergebnisse sind dann im "results" Ordner zu finden.
-Jede Simulation generiert eine "senkrechte.dat" und "waagrechte.dat" Datei für die Visualisierung mit 2dvis.py
     !!!Diese Dateien werden von meheren separaten Simulationen überschrieben und sind nicht kompatibel mit parallelen runs!!!
-Eine Datei entählt in den ersten 6 Zeilen die Simulationsparameter danach den Header und schließlich die Daten N N_h N_v
-Alle anderen Observablen und Mittelwerte wurden im Python file berechnet

--------------------------------------------
3. Kompilieren über Visual Studio
--------------------------------------------
-es wird Visual Studio benötigt
-Es muss die Datei "GCMC_VS.sln" (Projektmappe) in Visual Studio geöffnet werden welches die Datei beim starten automatisch kompiliert
-Beim starten in können die Parameter beim rechtsklicken von GCMC_VS im Projektmappen-Explorer 
 unter Eigenschaften -> Konfigurationseigenschaften -> Debugging -> Befehlsargumente angegeben werden:    
 <output_path + output_filename> <z> <thermalisation_steps> <interval> <total_measurements> <seed>
