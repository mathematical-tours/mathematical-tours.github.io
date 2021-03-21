function Faerbung = GraphColoringJohnson(Graph) 
%% Beschreibung:
% Diese Funktion liefert eine g?ltige Knotenf?rbung eines
% Graphen nach dem Johnson Algorithmus.
%
% Eingangswert: Ein Matlab Graph-Objekt. 
% Ausgabewert:  Ein Vektor der die F?rbung als Zahlenwerte
%               beinhaltet.
%
% Die Indizies des Ergebnis-Vektors entsprechen den Nummern
% der jeweiligen Knoten des Graph-Objektes. Die Zeit f?r die
% Berechnung der F?rbung wird mit ausgegeben.
%
% Beispiel:
%
% s = [1 1 2 2 3];
% t = [2 4 3 4 4];
% G = graph(s,t)
%
% G = 
%
%  graph with properties:
%
%    Edges: [5?1 table]
%    Nodes: [4?0 table]
%
%
% Farbvektor = GraphColoringJohnson(G)
% Elapsed time is 0.035626 seconds.
%
% Farbvektor =
%
%     0
%     1
%     0
%     2 

%% Initialisierung der Startwerte und Variablen

    % Zeiterfassung Start
    tic

    % Initialisierung der Startfarbe = 0
    Farbe_aktuell = -1;

    % Initialisiere die Faerbetabelle
    Faerbung = zeros(1, numnodes(Graph))';
    Faerbung(:) = Inf;

    % Knotennummerierung als Bezeichnung einf?hren
    Knotenbezeichnung = cell(numnodes(Graph),1);
    for i = 1:numnodes(Graph)
        Knotenbezeichnung{i} = num2str(i);
    end
    Graph.Nodes.Name = Knotenbezeichnung;

    % Definiere den Teilgraphen aus dem weiter gemacht wird
    Graph_neu = Graph;

    %%  Hier beginnt die Schleife

    % Wenn der aktuell betrachtete Untergraph abgebaut ist, erh?he die aktuelle
    % Farbe um 1, % baue den neuen Untergraaphen auf und f?hre den Johnson
    % Algorithmus solange aus, bis alle Knoten eine Farbe haben.
    while sum(isinf(Faerbung)) > 0

        % N?chste Farbe w?hlen
        Farbe_aktuell = Farbe_aktuell + 1;

        % Der neue Untergraph besteht aus allen Knoten des Original-Graphen
        % abz?glich der bereits gef?rbten Knoten. Diese Information wird aus
        % der akteullen Faerbungs-Liste gezogen
        Graph_neu = Graph;
        Graph_neu = rmnode(Graph_neu, find(~isinf(Faerbung)));

        % F?hre Johnson Algorithmus aus bis der aktuell betrachtet Untergraph
        % vollst?ndig abgebaut ist
        while numnodes(Graph_neu) > 0

            %W?hle den Knoten mit dem geringsten Grad
            [Grad_gewaehlter_Knoten, Index_Knoten_gewaehlt] = min(degree(Graph_neu, Graph_neu.Nodes.Name));

            % Finde die Nachbarn des Knotes mit dem geringsten Grad
            Nachbarn_akt_Knoten = neighbors(Graph_neu, (cell2mat(Graph_neu.Nodes.Name(Index_Knoten_gewaehlt))));

            % Setze f?r den gew?hlten Knoten die aktuell g?ltige Farbe
            Faerbung(str2double(cell2mat(Graph_neu.Nodes.Name(Index_Knoten_gewaehlt)))) = Farbe_aktuell;

            % Entferne den aktuell gew?hlten Knoten
            Graph_neu = rmnode(Graph_neu, Graph_neu.Nodes.Name(Index_Knoten_gewaehlt));

            % Entferne die Nachbarn des gef?rbten Knotens
            Graph_neu = rmnode(Graph_neu, Nachbarn_akt_Knoten);

        end

    end

    % Zeiterfassung Stop
    % toc

end