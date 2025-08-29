
%
function [ ] = drawBestTour(currentSolution , graph, fitness)
figure(2);
for i = 1 : length(currentSolution)-1
    
    currentNode = currentSolution(i);
    nextNode =  currentSolution(i+1);
    
    x1 = graph.node(currentNode).x;
    y1 = graph.node(currentNode).y;
    
    x2 = graph.node(nextNode).x;
    y2 = graph.node(nextNode).y;
    
    X = [x1 , x2];
    Y = [y1, y2];
    
    plot (X, Y, '-r');
    text(x1+0.2, y1,num2str(currentSolution(i)));
    hold on;

end


title(['Total length: ',num2str(fitness)]);
box('on');
hold off;

return
end

