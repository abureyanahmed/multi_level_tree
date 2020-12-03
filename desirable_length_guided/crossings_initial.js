function randomNumber(min, max) {  
    min = Math.ceil(min); 
    max = Math.floor(max); 
    return Math.floor(Math.random() * (max - min + 1)) + min; 
}

function straight_or_circular_force(cnt, speed, type)
{
    if(type=="straight")
    {
      return 1;
    }
    else
    {
      //return Math.cos(speed*cnt*Math.PI/360)+1;
      return Math.cos(speed*cnt*Math.PI/360)*.5+.5;
    }
}

function edgeInArray(removed_edges, u, v)
{
  for(var i=0;i<2;i++)
  {
    if(removed_edges[i][0]==u && removed_edges[i][1]==v)return true;
    if(removed_edges[i][0]==v && removed_edges[i][1]==u)return true;
  }
  return false;
}

function removeCrossingUsingRotation(cntr_vertex, other_vertex, crossing_line){
    //console.log("center vertex", id_to_label[cntr_vertex.id], cntr_vertex.x, cntr_vertex.y, "other vertex", id_to_label[other_vertex.id], other_vertex.x, other_vertex.y);
    var rotation_dir = 1;
    /*var rotation_dir = randomNumber(0, 1);
    if(rotation_dir==0)rotation_dir=1;
    else rotation_dir=-1;*/
    let number_of_angles = 100;
    var theta=Math.atan((other_vertex.y-cntr_vertex.y)/(other_vertex.x-cntr_vertex.x));
    for(var i=0;i<number_of_angles;i++)
    {
      let tx = other_vertex.x-cntr_vertex.x;
      let ty = other_vertex.y-cntr_vertex.y;
      //console.log("tx, ty", tx, ty);
      let rx = tx*Math.cos(theta)-ty*Math.sin(theta);
      let ry = tx*Math.sin(theta)+ty*Math.cos(theta);
      tx = rx+cntr_vertex.x;
      ty = ry+cntr_vertex.y;
      //console.log("tx, ty", tx, ty);
      var line1 = [
        [cntr_vertex.x, cntr_vertex.y],
        [tx, ty]
      ];

      var line2 = [
        [crossing_line[0].x, crossing_line[0].y],
        [crossing_line[1].x, crossing_line[1].y]
      ];

      if(linesCross(line1, line2)==false)return [tx, ty];
      theta=theta+(rotation_dir)*(Math.PI/number_of_angles);
    }
}

function removeAllCrossingsUsingRotation(cntr_vertex, other_vertex, inputLink){
    //console.log("center vertex", id_to_label[cntr_vertex.id], cntr_vertex.x, cntr_vertex.y, "other vertex", id_to_label[other_vertex.id], other_vertex.x, other_vertex.y);
    var oldX = other_vertex.x;
    var oldY = other_vertex.y;
    var rotation_dir = 1;
    /*var rotation_dir = randomNumber(0, 1);
    if(rotation_dir==0)rotation_dir=1;
    else rotation_dir=-1;*/
    let number_of_angles = 100;
    if(my_edges.length>2000)
    {
      number_of_angles = 1000;
    }
    var theta=Math.atan((other_vertex.y-cntr_vertex.y)/(other_vertex.x-cntr_vertex.x));
    for(var i=0;i<number_of_angles;i++)
    {
      let tx = other_vertex.x-cntr_vertex.x;
      let ty = other_vertex.y-cntr_vertex.y;
      //console.log("tx, ty", tx, ty);
      let rx = tx*Math.cos(theta)-ty*Math.sin(theta);
      let ry = tx*Math.sin(theta)+ty*Math.cos(theta);
      tx = rx+cntr_vertex.x;
      ty = ry+cntr_vertex.y;
      //console.log("tx, ty", tx, ty);
      other_vertex.x = tx;
      other_vertex.y = ty;
      var crossings = linkCrossingsWithInputLink(inputLink);
      //console.log(crossings);
      if(crossings.length==0)return [tx, ty];
      theta=theta+(rotation_dir)*(Math.PI/number_of_angles);
    }
    //console.log("Haven't found any crossing free position");
    return [oldX, oldY];
}

function get_child_nodes(firstEdge, secondEdge){
    var comp1Size = null;
    var comp1 = null;
    var removed_edges = [];
    myGraph.removeEdge(firstEdge[0].id, firstEdge[1].id);
    var part1 = myGraph.bfs(firstEdge[0].id);
    var part2 = myGraph.bfs(firstEdge[1].id);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {   
        if((part1[i]==secondEdge[0].id)||(part1[i]==secondEdge[1].id))
        { 
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {   
        var t = firstEdge[0];
        firstEdge[0] = firstEdge[1];
        firstEdge[1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp1Size = part1.length;
    comp1 = part1;
    myGraph.addEdge(firstEdge[0].id, firstEdge[1].id);
    return [comp1, firstEdge];
}

function update_drawing(arr, other_vertex, oldX, oldY, comp)
{
      let nx = arr[0];
      let ny = arr[1];
      graph.graphData.nodes[other_vertex.id].x = oldX;
      graph.graphData.nodes[other_vertex.id].y = oldY;
      let tx = nx-oldX;
      let ty = ny-oldY;
      for(var i=0;i<comp.length;i++)
      {
        graph.graphData.nodes[comp[i]].fx += tx;
        graph.graphData.nodes[comp[i]].fy += ty;
      }
      d3.select('#nodes_testSvgId').selectAll('g').attr("transform", d => {return 'translate('+[d.x,d.y]+')';});
      d3.select('#links_testSvgId').selectAll('line').attr('x1', d => d.source.x).attr('y1', d => d.source.y).attr('x2', d => d.target.x).attr('y2', d => d.target.y).style("stroke", d => d.stroke).style("stroke-width", d => d["stroke-width"]);
}

function remove_existing_crossings(){
  var res = linkCrossingsParam(graph.graphData.links);
  var old_length = 0;
  while(res.length>0)
  {
    //var crossingIndex = randomNumber(0, res.length-1);
    for(var crossingIndex=0;crossingIndex<res.length;crossingIndex++){
      //console.log("crossingIndex", crossingIndex);
      var crossingPair = res[crossingIndex];
      var arr = get_child_nodes(crossingPair[0], crossingPair[1]);
      var comp1 = arr[0];
      crossingPair[0] = arr[1];
      //console.log(crossingPair[0]);
      var arr = get_child_nodes(crossingPair[1], crossingPair[0]);
      var comp2 = arr[0];
      crossingPair[1] = arr[1];
      //console.log(crossingPair[1]);
      var new_coord = null;
      if(comp1.length<comp2.length)
        {
          var oldX = crossingPair[0][0].x;
          var oldY = crossingPair[0][0].y;
          new_coord = removeAllCrossingsUsingRotation(crossingPair[0][1], crossingPair[0][0], graph.graphData.links[crossingPair[0][2]]);
          update_drawing(new_coord, crossingPair[0][0], oldX, oldY, comp1);
          //console.log("Center", crossingPair[0][1], "other", crossingPair[0][0]);
          if((oldX==new_coord[0])&&(oldX==new_coord[1]))
          {
            edge_distance[crossingPair[0][2]] = Math.max(10, edge_distance[crossingPair[0][2]]-5)
          }
        }
      else
        {
          var oldX = crossingPair[1][0].x;
          var oldY = crossingPair[1][0].y;
          new_coord = removeAllCrossingsUsingRotation(crossingPair[1][1], crossingPair[1][0], graph.graphData.links[crossingPair[1][2]]);
          update_drawing(new_coord, crossingPair[1][0], oldX, oldY, comp2);
          //console.log("Center", crossingPair[1][1], "other", crossingPair[1][0]);
          if((oldX==new_coord[0])&&(oldX==new_coord[1]))
          {
            edge_distance[crossingPair[0][2]] = Math.max(10, edge_distance[crossingPair[0][2]]-5)
          }
        }
    }
    res = linkCrossingsParam(graph.graphData.links);
    if(old_length==res.length){
      if(res.length>0){
        //addEdgePhase = false;
        //removeSmallCompIntervalMultipleSafe = setInterval(removeSmallCompAutoMultipleSafe, intervalTime);
        console.log('Unable to remove crossing, leaving as it is.');
      }
      break;
    }
    old_length = res.length;
    //document.getElementById("number_of_crossings").value = res.length;
    console.log("crossings", res.length);
  }
}

var removeSmallCompMultiple = function (crossingIndex, res) {
    // select a crossing pair randomly
    //var crossingIndex = randomNumber(0, res.length-1);
    var crossingPair = res[crossingIndex];

    var comp1Size = null;
    var comp1 = null;
    var removed_edges = [];
    myGraph.removeEdge(crossingPair[0][0].id, crossingPair[0][1].id);
    removed_edges.push([crossingPair[0][0].id, crossingPair[0][1].id]);
    var part1 = myGraph.bfs(crossingPair[0][0].id);
    var part2 = myGraph.bfs(crossingPair[0][1].id);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {
        if((part1[i]==crossingPair[1][0].id)||(part1[i]==crossingPair[1][1].id))
        {
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {
        var t = crossingPair[0][0];
        crossingPair[0][0] = crossingPair[0][1];
        crossingPair[0][1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp1Size = part1.length;
    comp1 = part1;
    var comp2Size = null;
    var comp2 = null;
    myGraph.removeEdge(crossingPair[1][0].id, crossingPair[1][1].id);
    removed_edges.push([crossingPair[1][0].id, crossingPair[1][1].id]);
    var part1 = myGraph.bfs(crossingPair[1][0].id);
    var part2 = myGraph.bfs(crossingPair[1][1].id);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {
        if((part1[i]==crossingPair[0][0].id)||(part1[i]==crossingPair[0][1].id))
        {
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {
        var t = crossingPair[1][0];
        crossingPair[1][0] = crossingPair[1][1];
        crossingPair[1][1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp2Size = part1.length;
    comp2 = part1;
    var commonCompSize = part2.length;
    var commonComp = part2;
    myGraph.addEdge(crossingPair[0][0].id, crossingPair[0][1].id);
    if(!edgeInArray(removed_edges, crossingPair[0][0].id, crossingPair[0][1].id))console.log(crossingPair[0][0].id, crossingPair[0][1].id, "not in", removed_edges);
    myGraph.addEdge(crossingPair[1][0].id, crossingPair[1][1].id);
    if(!edgeInArray(removed_edges, crossingPair[1][0].id, crossingPair[1][1].id))console.log(crossingPair[1][0].id, crossingPair[1][1].id, "not in", removed_edges);
    // In both edge the first vertex contains the separate component, place force accordingly
    //console.log("comp1Size:", comp1Size, "comp2Size:", comp2Size, "commonCompSize:", commonCompSize);

    let force_factor = straight_or_circular_force(crossing_step_cnt, crossing_step_speed, crossing_step_type);
    crossing_step_cnt = crossing_step_cnt + 1;
    /*if((comp1Size==1) || (comp2Size==1))
    {
      let cntr_vertex = -1;
      let other_vertex = -1;
      let crossing_line = -1;
      if(comp1Size==1)
      {
        //console.log(id_to_label[comp1[0]], "is empty");
        if(crossingPair[0][0].id==comp1[0])
        {
          other_vertex=crossingPair[0][0];
          cntr_vertex=crossingPair[0][1];
        }
        else
        {
          other_vertex=crossingPair[0][1];
          cntr_vertex=crossingPair[0][0];
        }
        crossing_line = crossingPair[1];
      }
      else if(comp2Size==1)
      {
        //console.log(id_to_label[comp2[0]], "is empty");
        if(crossingPair[1][0].id==comp2[0])
        {
          other_vertex=crossingPair[1][0];
          cntr_vertex=crossingPair[1][1];
        }
        else
        {
          other_vertex=crossingPair[1][1];
          cntr_vertex=crossingPair[1][0];
        }
        crossing_line = crossingPair[0];
      }
      let arr = removeCrossingUsingRotation(cntr_vertex, other_vertex, crossing_line);
      let nx = arr[0];
      let ny = arr[1];
      //console.log(graph.graphData.nodes[other_vertex.id].x);
      graph.graphData.nodes[other_vertex.id].x = nx;
      graph.graphData.nodes[other_vertex.id].y = ny;
      //console.log(graph.graphData.nodes[other_vertex.id].x);
      d3.select('#nodes_testSvgId').selectAll('g').attr("transform", d => {return 'translate('+[d.x,d.y]+')';});
      d3.select('#links_testSvgId').selectAll('line').attr('x1', d => d.source.x).attr('y1', d => d.source.y).attr('x2', d => d.target.x).attr('y2', d => d.target.y).style("stroke", d => d.stroke).style("stroke-width", d => d["stroke-width"]);
    }
    else*/ 
    if((comp1Size<=comp2Size) && (comp1Size<=commonCompSize))
    {
      var moveAmountUp = crossingPair[0][1].x-crossingPair[0][0].x;
      var moveAmountRight = crossingPair[0][1].y-crossingPair[0][0].y;
      var nodeId = crossingPair[0][1].id;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight*force_factor, moveAmountUp*force_factor);
      moveParMultiple(new Set(comp1), moveAmountUp*force_factor, moveAmountRight*force_factor);
    }
    else if((comp2Size<=comp1Size) && (comp2Size<=commonCompSize))
    {
      var moveAmountUp = crossingPair[1][1].x-crossingPair[1][0].x;
      var moveAmountRight = crossingPair[1][1].y-crossingPair[1][0].y;
      var nodeId = crossingPair[1][1].id;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight*force_factor, moveAmountUp*force_factor);
      moveParMultiple(new Set(comp2), moveAmountUp*force_factor, moveAmountRight*force_factor);
    }
    else
    {
      var moveAmountUp = crossingPair[0][0].x-crossingPair[0][1].x;
      var moveAmountRight = crossingPair[0][0].y-crossingPair[0][1].y;
      var nodeId = crossingPair[0][0].id;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight*force_factor, moveAmountUp*force_factor);
      moveParMultiple(new Set(comp1), moveAmountUp*force_factor, moveAmountRight*force_factor);
    }

}

var removeSmallCompMultipleOld = function (crossingIndex, res) {
    // select a crossing pair randomly
    //var crossingIndex = randomNumber(0, res.length-1);
    var crossingPair = res[crossingIndex];
    console.log(crossingPair);

    var comp1Size = null;
    var comp1 = null;
    //myGraph.removeEdge(crossingPair[0][0].index, crossingPair[0][1].index);
    myGraph.removeEdge(crossingPair[0][0].id, crossingPair[0][1].id);
    console.log("removing edge:", crossingPair[0][0].index, crossingPair[0][1].index);
    console.log("Adj list:", myGraph.AdjList);
    var part1 = myGraph.bfs(crossingPair[0][0].index);
    var part2 = myGraph.bfs(crossingPair[0][1].index);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {
        if((part1[i]==crossingPair[1][0].index)||(part1[i]==crossingPair[1][1].index))
        {
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {
        var t = crossingPair[0][0];
        crossingPair[0][0] = crossingPair[0][1];
        crossingPair[0][1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp1Size = part1.length;
    comp1 = part1;
    var comp2Size = null;
    var comp2 = null;
    //myGraph.removeEdge(crossingPair[1][0].index, crossingPair[1][1].index);
    myGraph.removeEdge(crossingPair[1][0].id, crossingPair[1][1].id);
    console.log("removing edge:", crossingPair[1][0].index, crossingPair[1][1].index);
    console.log("Adj list:", myGraph.AdjList);
    var part1 = myGraph.bfs(crossingPair[1][0].index);
    var part2 = myGraph.bfs(crossingPair[1][1].index);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {
        if((part1[i]==crossingPair[0][0].index)||(part1[i]==crossingPair[0][1].index))
        {
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {
        var t = crossingPair[1][0];
        crossingPair[1][0] = crossingPair[1][1];
        crossingPair[1][1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp2Size = part1.length;
    comp2 = part1;
    var commonCompSize = part2.length;
    var commonComp = part2;
    //myGraph.addEdge(crossingPair[0][0].index, crossingPair[0][1].index);
    myGraph.addEdge(crossingPair[0][0].id, crossingPair[0][1].id);
    console.log("adding edge:", crossingPair[0][0].index, crossingPair[0][1].index);
    //myGraph.addEdge(crossingPair[1][0].index, crossingPair[1][1].index);
    myGraph.addEdge(crossingPair[1][0].id, crossingPair[1][1].id);
    console.log("adding edge:", crossingPair[1][0].index, crossingPair[1][1].index);
    // In both edge the first vertex contains the separate component, place force accordingly
    //console.log("comp1Size:", comp1Size, "comp2Size:", comp2Size, "commonCompSize:", commonCompSize);

    let force_factor = straight_or_circular_force(crossing_step_cnt, crossing_step_speed, crossing_step_type);
    crossing_step_cnt = crossing_step_cnt + 1;
    if((comp1Size<=comp2Size) && (comp1Size<=commonCompSize))
    {
      var moveAmountUp = crossingPair[0][1].x-crossingPair[0][0].x;
      var moveAmountRight = crossingPair[0][1].y-crossingPair[0][0].y;
      var nodeId = crossingPair[0][1].index;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight*force_factor, moveAmountUp*force_factor);
      moveParMultiple(new Set(comp1), moveAmountUp*force_factor, moveAmountRight*force_factor);
    }
    else if((comp2Size<=comp1Size) && (comp2Size<=commonCompSize))
    {
      var moveAmountUp = crossingPair[1][1].x-crossingPair[1][0].x;
      var moveAmountRight = crossingPair[1][1].y-crossingPair[1][0].y;
      var nodeId = crossingPair[1][1].index;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight*force_factor, moveAmountUp*force_factor);
      moveParMultiple(new Set(comp2), moveAmountUp*force_factor, moveAmountRight*force_factor);
    }
    else
    {
      var moveAmountUp = crossingPair[0][0].x-crossingPair[0][1].x;
      var moveAmountRight = crossingPair[0][0].y-crossingPair[0][1].y;
      var nodeId = crossingPair[0][0].index;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight*force_factor, moveAmountUp*force_factor);
      moveParMultiple(new Set(comp1), moveAmountUp*force_factor, moveAmountRight*force_factor);
    }

}

var removeSmallCompAutoMultiple = function () {
  var res = linkCrossingsParam(graph.graphData.links);
  var crossingIndex = randomNumber(0, res.length-1);
  //var crossingIndex = 0;
  if(res.length>0){
    removeSmallCompMultiple(crossingIndex, res);
  }
  else
  {
    stopRemovingSmallCompMultiple();
    var t1 = performance.now();
    console.log("Call to doSomething took " + (t1 - t0) + " milliseconds.");
  }
}

function stopRemovingSmallCompMultiple() {
  clearInterval(removeSmallCompIntervalMultiple);
}

function unhighlightEdge(){
    if(previousCrossingIndex!=null)
    {
      graph.graphData.links[previousCrossingIndex].stroke = "#76963e";
      graph.graphData.links[previousCrossingIndex]["stroke-width"] = "1.0px";
    }
}

function moveAwayFromCenter(forceFactor){
  var nodeDepth = myGraph.bfs_depth(0);
  //console.log(nodeDepth);
  var i=0;
  for(var i=1;i<nodeDepth.length;i++)
  {
    var vx = graph.graphData.nodes[i].x - graph.graphData.nodes[0].x;
    var vy = graph.graphData.nodes[i].y - graph.graphData.nodes[0].y;
    var nrm = Math.sqrt(vx*vx + vy*vy);
    vx = vx/nrm;
    vy = vy/nrm;
    movePar(i, vx*forceFactor, vy*forceFactor);
  }
  //console.log("Added force from the direction from center");
}

previousCrossingIndex = null;
var removeSmallCompAutoMultipleSafe = function () {
  moveAwayFromCenter(100);
  var res = linkCrossingsParam(graph.graphData.links);
  //var crossingIndex = randomNumber(0, res.length-1);
  var crossingIndex = 0;
  document.getElementById("number_of_crossings").value = res.length;
  //if(res.length>0){
  //if(res.length>50){
  if(!addEdgePhase){
    let edgeIndex = res[crossingIndex][0][2];
    graph.graphData.links[edgeIndex].stroke = "red";
    graph.graphData.links[edgeIndex]["stroke-width"] = "5.0px";
    if(edgeIndex!=previousCrossingIndex)unhighlightEdge();
    previousCrossingIndex = edgeIndex;
    removeSmallCompMultiple(crossingIndex, res);
  }
  else
  {
    stopRemovingSmallCompMultipleSafe();
    unhighlightEdge();
    addEdgeInterval = setInterval(startAddingEdges, intervalTime);
  }
}

function stopRemovingSmallCompMultipleSafe() {
  clearInterval(removeSmallCompIntervalMultipleSafe);
}

var removeSmallComp = function (crossingIndex, res, forceFactor) {
    // select a crossing pair randomly
    //var crossingIndex = randomNumber(0, res.length-1);
    var crossingPair = res[crossingIndex];

    var comp1Size = null;
    myGraph.removeEdge(crossingPair[0][0].index, crossingPair[0][1].index);
    var part1 = myGraph.bfs(crossingPair[0][0].index);
    var part2 = myGraph.bfs(crossingPair[0][1].index);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {
        if((part1[i]==crossingPair[1][0].index)||(part1[i]==crossingPair[1][1].index))
        {
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {
        var t = crossingPair[0][0];
        crossingPair[0][0] = crossingPair[0][1];
        crossingPair[0][1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp1Size = part1.length;
    var comp2Size = null;
    myGraph.removeEdge(crossingPair[1][0].index, crossingPair[1][1].index);
    var part1 = myGraph.bfs(crossingPair[1][0].index);
    var part2 = myGraph.bfs(crossingPair[1][1].index);
    var part1ContainsSecondEdge = false;
    for(var i=0; i<part1.length; i++)
    {
        if((part1[i]==crossingPair[0][0].index)||(part1[i]==crossingPair[0][1].index))
        {
          part1ContainsSecondEdge = true;
          break;
        }
    }
    if(part1ContainsSecondEdge)
    {
        var t = crossingPair[1][0];
        crossingPair[1][0] = crossingPair[1][1];
        crossingPair[1][1] = t;
        t = part1;
        part1 = part2;
        part2 = t;
    }
    comp2Size = part1.length;
    var commonCompSize = part2.length;
    myGraph.addEdge(crossingPair[0][0].index, crossingPair[0][1].index);
    myGraph.addEdge(crossingPair[1][0].index, crossingPair[1][1].index);
    // In both edge the first vertex contains the separate component, place force accordingly
    //console.log("comp1Size:", comp1Size, "comp2Size:", comp2Size, "commonCompSize:", commonCompSize);

    if((comp1Size<=comp2Size) && (comp1Size<=commonCompSize))
    {
      var moveAmountUp = crossingPair[0][1].x-crossingPair[0][0].x;
      var moveAmountRight = crossingPair[0][1].y-crossingPair[0][0].y;
      var nodeId = crossingPair[0][1].index;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight, moveAmountUp);
      movePar(nodeId, moveAmountUp*forceFactor, moveAmountRight*forceFactor);
    }
    else if((comp2Size<=comp1Size) && (comp2Size<=commonCompSize))
    {
      var moveAmountUp = crossingPair[1][1].x-crossingPair[1][0].x;
      var moveAmountRight = crossingPair[1][1].y-crossingPair[1][0].y;
      var nodeId = crossingPair[1][1].index;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight, moveAmountUp);
      movePar(nodeId, moveAmountUp*forceFactor, moveAmountRight*forceFactor);
    }
    else
    {
      var moveAmountUp = crossingPair[0][0].x-crossingPair[0][1].x;
      var moveAmountRight = crossingPair[0][0].y-crossingPair[0][1].y;
      var nodeId = crossingPair[0][0].index;
      //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight, moveAmountUp);
      movePar(nodeId, moveAmountUp*forceFactor, moveAmountRight*forceFactor);
    }

}

prevResLength = null;
forceFactor = 1;
var removeSmallCompAuto = function () {
  var res = linkCrossingsParam(graph.graphData.links);
  var crossingIndex = randomNumber(0, res.length-1);
  //var crossingIndex = 0;
  if(prevResLength!=null){
    //if(prevResLength>=res.length){if(forceFactor<128)forceFactor *=2;}
    if(prevResLength>=res.length){if(forceFactor<32)forceFactor *=2;}
    else{forceFactor = 1;}
  }
  else{prevResLength = res.length;}
  if(res.length>0){
    removeSmallComp(crossingIndex, res, forceFactor);
  }
  else
  {
    stopRemovingSmallComp();
    var t1 = performance.now();
    console.log("Call to doSomething took " + (t1 - t0) + " milliseconds.");
  }
}

function stopRemovingSmallComp() {
  clearInterval(removeSmallCompInterval);
}

//var initialize = function () {
var addCrossingForce = function () {
  //console.log(graph.graphData.links);
  links = graph.graphData.links;
  m = links.length;
  res = linkCrossings();
  //console.log("Number of crossings:", res.length);
  if(res.length>0){
    // select a crossing pair randomly
    var crossingIndex = randomNumber(0, res.length-1);
    var crossingPair = res[crossingIndex];
    // swap the edges with probability .5
    var coin = randomNumber(0, 1);
    if(coin==0){
      var t = crossingPair[0];
      crossingPair[0] = crossingPair[1];
      crossingPair[1] = t;
    }
    // swap the points with probability .5
    var coin = randomNumber(0, 1);
    if(coin==0){
      var t = crossingPair[0][0];
      crossingPair[0][0] = crossingPair[0][1];
      crossingPair[0][1] = t;
    }
    var moveAmountUp = crossingPair[0][0].x-graph.graphData.nodes[0].x;
    var moveAmountRight = crossingPair[0][0].y-graph.graphData.nodes[0].y;
    var nodeId = crossingPair[0][1].index;
    //console.log("Moving ", id_to_label[nodeId], " to ", moveAmountRight, moveAmountUp);
    movePar(nodeId, moveAmountUp, moveAmountRight);
  }
  else{
    stopCrossingForce();
    var t1 = performance.now();
    console.log("Call to doSomething took " + (t1 - t0) + " milliseconds.");
  }
}

function stopCrossingForce() {
  clearInterval(crossingForceInterval);
}

function direction (pi, pj, pk) {
    var p1 = [pk[0] - pi[0], pk[1] - pi[1]];
    var p2 = [pj[0] - pi[0], pj[1] - pi[1]];
    return p1[0] * p2[1] - p2[0] * p1[1];
}

// Is point k on the line segment formed by points i and j?
// Inclusive, so if pk == pi or pk == pj then return true.
function onSegment (pi, pj, pk) {
    return Math.min(pi[0], pj[0]) <= pk[0] &&
      pk[0] <= Math.max(pi[0], pj[0]) &&
      Math.min(pi[1], pj[1]) <= pk[1] &&
      pk[1] <= Math.max(pi[1], pj[1]);
}

function linesCross (line1, line2) {
    var d1, d2, d3, d4;

    // CLRS 2nd ed. pg. 937
    d1 = direction(line2[0], line2[1], line1[0]);
    d2 = direction(line2[0], line2[1], line1[1]);
    d3 = direction(line1[0], line1[1], line2[0]);
    d4 = direction(line1[0], line1[1], line2[1]);

    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
      ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
      return true;
    } else if (d1 === 0 && onSegment(line2[0], line2[1], line1[0])) {
      return true;
    } else if (d2 === 0 && onSegment(line2[0], line2[1], line1[1])) {
      return true;
    } else if (d3 === 0 && onSegment(line1[0], line1[1], line2[0])) {
      return true;
    } else if (d4 === 0 && onSegment(line1[0], line1[1], line2[1])) {
      return true;
    }

    return false;
}

function linksCross (link1, link2) {
    // Self loops are not intersections
    if (link1.index === link2.index ||
      link1.source === link1.target ||
      link2.source === link2.target) {
      return false;
    }

    // Links cannot intersect if they share a node
    if (link1.source === link2.source ||
      link1.source === link2.target ||
      link1.target === link2.source ||
      link1.target === link2.target) {
      return false;
    }

    var line1 = [
      [link1.source.x, link1.source.y],
      [link1.target.x, link1.target.y]
    ];

    var line2 = [
      [link2.source.x, link2.source.y],
      [link2.target.x, link2.target.y]
    ];

    return linesCross(line1, line2);
}

function linksCrossWithCrds (link1, link2, crd_x, crd_y) {
    // Self loops are not intersections
    if (link1.index === link2.index ||
      link1.source === link1.target ||
      link2.source === link2.target) {
      return false;
    }

    // Links cannot intersect if they share a node
    if (link1.source === link2.source ||
      link1.source === link2.target ||
      link1.target === link2.source ||
      link1.target === link2.target) {
      return false;
    }

    var line1 = [
      [crd_x[link1.source.id], crd_y[link1.source.id]],
      [crd_x[link1.target.id], crd_y[link1.target.id]]
    ];

    var line2 = [
      [crd_x[link2.source.id], crd_y[link2.source.id]],
      [crd_x[link2.target.id], crd_y[link2.target.id]]
    ];

    return linesCross(line1, line2);
}

function linkCrossings () {
    var i, j, c = 0, link1, link2, line1, line2;;
    var res = [];

    // Sum the upper diagonal of the edge crossing matrix.
    for (i = 0; i < m; ++i) {
      for (j = i + 1; j < m; ++j) {
        link1 = links[i], link2 = links[j];

        // Check if link i and link j intersect
        if (linksCross(link1, link2)) {
          line1 = [
            [link1.source.x, link1.source.y],
            [link1.target.x, link1.target.y]
          ];
          line2 = [
            [link2.source.x, link2.source.y],
            [link2.target.x, link2.target.y]
          ];
          ++c;
          //console.log(id_to_label[link1.source.index], ",", id_to_label[link1.target.index], " crosses ", id_to_label[link2.source.index], ",", id_to_label[link2.target.index]);
          res.push([[link1.source, link1.target], [link2.source, link2.target]]);
          //d += Math.abs(idealAngle - acuteLinesAngle(line1, line2));
        }
      }
    }

    return res;
}

function linkCrossingsWithCrds (crd_x, crd_y) {
    var i, j, c = 0, link1, link2, line1, line2;;
    var res = [];

    var links = graph.graphData.links;
    let m = links.length;
    // Sum the upper diagonal of the edge crossing matrix.
    for (i = 0; i < m; ++i) {
      for (j = i + 1; j < m; ++j) {
        link1 = links[i], link2 = links[j];

        // Check if link i and link j intersect
        if (linksCrossWithCrds(link1, link2, crd_x, crd_y)) {
          /*line1 = [
            [link1.source.x, link1.source.y],
            [link1.target.x, link1.target.y]
          ];
          line2 = [
            [link2.source.x, link2.source.y],
            [link2.target.x, link2.target.y]
          ];
          ++c;*/
          //console.log(id_to_label[link1.source.index], ",", id_to_label[link1.target.index], " crosses ", id_to_label[link2.source.index], ",", id_to_label[link2.target.index]);
          res.push([[link1.source, link1.target], [link2.source, link2.target]]);
          //d += Math.abs(idealAngle - acuteLinesAngle(line1, line2));
        }
      }
    }

    return res;
}

function linkCrossingsWithInputLink (inputLink) {
    var i, j, c = 0, link1, link2, line1, line2;;
    var res = [];

    // Sum the upper diagonal of the edge crossing matrix.
    var links = graph.graphData.links;
    let m = links.length;
    for (var i = 0; i < m; ++i) {
        link1 = links[i], link2 = inputLink;

        // Check if link i and link j intersect
        if (linksCross(link1, link2)) {
          line1 = [
            [link1.source.x, link1.source.y],
            [link1.target.x, link1.target.y]
          ];
          line2 = [
            [link2.source.x, link2.source.y],
            [link2.target.x, link2.target.y]
          ];
          ++c;
          //console.log(id_to_label[link1.source.index], ",", id_to_label[link1.target.index], " crosses ", id_to_label[link2.source.index], ",", id_to_label[link2.target.index]);
          res.push([[link1.source, link1.target], [link2.source, link2.target]]);
          //d += Math.abs(idealAngle - acuteLinesAngle(line1, line2));
        }
    }

    return res;
}

function hasLinkCrossingsWithInputLink (inputLink, crd_x, crd_y) {
    var i, j, c = 0, link1, link2, line1, line2;;
    //var res = [];

    // Sum the upper diagonal of the edge crossing matrix.
    var links = graph.graphData.links;
    let m = links.length;
    for (var i = 0; i < m; ++i) {
        link1 = links[i], link2 = inputLink;

        // Check if link i and link j intersect
        if (linksCrossWithCrds(link1, link2, crd_x, crd_y)) {
          /*line1 = [
            [link1.source.x, link1.source.y],
            [link1.target.x, link1.target.y]
          ];
          line2 = [
            [link2.source.x, link2.source.y],
            [link2.target.x, link2.target.y]
          ];
          ++c;*/
          //console.log(id_to_label[link1.source.index], ",", id_to_label[link1.target.index], " crosses ", id_to_label[link2.source.index], ",", id_to_label[link2.target.index]);
          //res.push([[link1.source, link1.target], [link2.source, link2.target]]);
          //d += Math.abs(idealAngle - acuteLinesAngle(line1, line2));
          return true;
        }
    }

    //return res;
    return false;
}

function linkCrossingsParam (links) {
    var i, j, c = 0, link1, link2, line1, line2;;
    var res = [];

    // Sum the upper diagonal of the edge crossing matrix.
    let m = links.length;
    for (i = 0; i < m; ++i) {
      for (j = i + 1; j < m; ++j) {
        link1 = links[i], link2 = links[j];

        // Check if link i and link j intersect
        if (linksCross(link1, link2)) {
          line1 = [
            [link1.source.x, link1.source.y],
            [link1.target.x, link1.target.y]
          ];
          line2 = [
            [link2.source.x, link2.source.y],
            [link2.target.x, link2.target.y]
          ];
          ++c;
          //console.log(id_to_label[link1.source.index], ",", id_to_label[link1.target.index], " crosses ", id_to_label[link2.source.index], ",", id_to_label[link2.target.index]);
          res.push([[link1.source, link1.target, link1.index], [link2.source, link2.target, link2.index]]);
          //d += Math.abs(idealAngle - acuteLinesAngle(line1, line2));
        }
      }
    }

    return res;
}

//overall the idea does not work
//reason one: there is no efficient way to calculate all points in a region
//reason two, even you detect all the points in a region containing the edge, there might be a vertex far away from that region but still crosses the edge
function search(quadtree, xmin, ymin, xmax, ymax) {
  const results = [];
  quadtree.visit(function(node, x1, y1, x2, y2) {
    if (!node.length) {
      do {
        var d = node.data;
        //if (d[0] >= xmin && d[0] < xmax && d[1] >= ymin && d[1] < ymax) {
        //if (d.x >= xmin && d.x < xmax && d.y >= ymin && d.y < ymax) {
        if (crd_x[d.id] >= xmin && crd_x[d.id] < xmax && crd_y[d.id] >= ymin && crd_y[d.id] < ymax) {
          results.push(d);
        }
      } while (node = node.next);
    }
    return x1 >= xmax || y1 >= ymax || x2 < xmin || y2 < ymin;
  });
  return results;
}

function linkCrossingsParamQuadTree (link1, crd_x, crd_y, node_to_links, quadtree, max_edge_len) {
  /*let links = graph.graphData.links;
  var quadtree = d3.quadtree()
  .x(function (d) { return d.x; })
  .y(function (d) { return d.y; });
  var data = [links[0].source];
  for(var i=0;i<links.length;i++)
    data.push(links[i].target);
  quadtree.addAll(data);
  let max_edge_len = 500;
  for(var i = 0;i<links.length;i++)
  {
    let dis = Math.sqrt((links[i].source.x-links[i].target.x)*(links[i].source.x-links[i].target.x) + (links[i].source.y-links[i].target.y)*(links[i].source.y-links[i].target.y));
    if(max_edge_len<dis)max_edge_len=dis;
  }
  let count = 0;
  for(var i = 0;i<links.length;i++)
  {
    let xmin = Math.min(links[i].source.x, links[i].target.x)-max_edge_len;
    let ymin = Math.min(links[i].source.y, links[i].target.y)-max_edge_len;
    let xmax = Math.max(links[i].source.x, links[i].target.x)+max_edge_len;
    let ymax = Math.max(links[i].source.y, links[i].target.y)+max_edge_len;
    let crossings = search(quadtree, xmin, ymin, xmax, ymax);
    count = count + crossings.length;
  }*/
  let xmin = Math.min(crd_x[link1.source.id], crd_x[link1.target.id])-max_edge_len;
  let ymin = Math.min(crd_y[link1.source.id], crd_y[link1.target.id])-max_edge_len;
  let xmax = Math.max(crd_x[link1.source.id], crd_x[link1.target.id])+max_edge_len;
  let ymax = Math.max(crd_y[link1.source.id], crd_y[link1.target.id])+max_edge_len;
  let nbrNodes = search(quadtree, xmin, ymin, xmax, ymax);
  for(var i=0;i<nbrNodes.length;i++)
  {
     let id = nbrNodes[i].id;
     for(let j=0;j<node_to_links[id].length;j++)
     {
       let link2 = node_to_links[id][j];
       if (linksCrossWithCrds(link1, link2, crd_x, crd_y))
       {
         return true;
       }
     }
  }
  return false;
}


