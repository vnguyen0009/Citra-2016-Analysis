<?xml version="1.0" encoding="UTF-8"?><project xmlns="http://grogra.de/registry" graph="graph.xml">
 <import plugin="de.grogra.imp" version="1.5"/>
 <import plugin="de.grogra.pf" version="1.5"/>
 <import plugin="de.grogra.imp3d" version="1.5"/>
 <import plugin="de.grogra.rgg" version="1.5"/>
 <import plugin="de.grogra.vecmath" version="1.5"/>
 <import plugin="de.grogra.math" version="1.5"/>
 <registry>
  <ref name="project">
   <ref name="objects">
    <ref name="files">
     <de.grogra.pf.ui.registry.SourceFile mimeType="text/x-grogra-rgg" name="pfs:NAR.rgg"/>
    </ref>
    <ref name="datasets">
     <de.grogra.pf.registry.SharedValue name="dat" type="de.grogra.pf.data.Dataset" value="serialized:rO0ABXNyABlkZS5ncm9ncmEucGYuZGF0YS5EYXRhc2V0Pri86I/xbusCAAlaAAxzZXJpZXNJblJvd3NMAARiaW5zdAAVTGphdmEvdXRpbC9BcnJheUxpc3Q7TAANY2F0ZWdvcnlMYWJlbHQAEkxqYXZhL2xhbmcvU3RyaW5nO0wACmNvbHVtbktleXNxAH4AAUwABW9yZGVydAAcTG9yZy9qZnJlZS9kYXRhL0RvbWFpbk9yZGVyO0wAB3Jvd0tleXNxAH4AAUwABHJvd3NxAH4AAUwABXRpdGxlcQB+AAJMAAp2YWx1ZUxhYmVscQB+AAJ4cABzcgATamF2YS51dGlsLkFycmF5TGlzdHiB0h2Zx2GdAwABSQAEc2l6ZXhwAAAAAHcEAAAAAHhwc3EAfgAFAAAABHcEAAAABHQADmN1cnJlbnROb2RlTnVtdAAJQnJhbmNoTnVtdAAIVE9UQUxMRU50AAlNU0xFQUZOVU14cHNxAH4ABQAAAAB3BAAAAAB4c3EAfgAFAAAAAHcEAAAAAHh0AANkYXRw"/>
    </ref>
    <ref name="meta">
     <de.grogra.pf.registry.NodeReference name="NAR" ref="23794"/>
    </ref>
   </ref>
  </ref>
  <ref name="workbench">
   <ref name="state">
    <de.grogra.pf.ui.registry.Layout name="layout">
     <de.grogra.pf.ui.registry.MainWindow>
      <de.grogra.pf.ui.registry.Split location="0.47713864">
       <de.grogra.pf.ui.registry.Split location="0.79613096" orientation="0">
        <de.grogra.pf.ui.registry.Split orientation="0">
         <de.grogra.pf.registry.Link source="/ui/panels/rgg/toolbar"/>
         <de.grogra.pf.ui.registry.PanelFactory source="/ui/panels/3d/defaultview">
          <de.grogra.pf.registry.Option name="panelId" type="java.lang.String" value="/ui/panels/3d/defaultview"/>
          <de.grogra.pf.registry.Option name="panelTitle" type="java.lang.String" value="View"/>
          <de.grogra.pf.registry.Option name="view" type="de.grogra.imp3d.View3D" value="graphDescriptor=[de.grogra.imp.ProjectGraphDescriptor]visibleScales={true true true true true true true true true true true true true true true}visibleLayers={true true true true true true true true true true true true true true true true}epsilon=1.0E-6 visualEpsilon=0.01 magnitude=1.0 camera=(minZ=0.1 maxZ=2000.0 projection=[de.grogra.imp3d.PerspectiveProjection aspect=1.0 fieldOfView=0.26838842]transformation=(-0.904451365790017 -0.4265767538445758 0.0 0.195950174505242 0.06391476844177667 -0.1355155879693838 0.9887114987661146 -20.325853981319593 -0.421761341632451 0.8942414654313142 0.149831813069368 -400.86086970717395 0.0 0.0 0.0 1.0))navigator=null"/>
         </de.grogra.pf.ui.registry.PanelFactory>
        </de.grogra.pf.ui.registry.Split>
        <de.grogra.pf.ui.registry.Split orientation="0">
         <de.grogra.pf.ui.registry.Tab selectedIndex="0">
          <de.grogra.pf.registry.Link source="/ui/panels/fileexplorer"/>
          <de.grogra.pf.registry.Link source="/ui/panels/objects/meta"/>
         </de.grogra.pf.ui.registry.Tab>
         <de.grogra.pf.registry.Link source="/ui/panels/statusbar"/>
        </de.grogra.pf.ui.registry.Split>
       </de.grogra.pf.ui.registry.Split>
       <de.grogra.pf.ui.registry.Split location="0.57738096" orientation="0">
        <de.grogra.pf.ui.registry.Tab selectedIndex="0">
         <de.grogra.pf.ui.registry.PanelFactory source="/ui/panels/texteditor">
          <de.grogra.pf.registry.Option name="documents" type="java.lang.String" value="&quot;\&quot;pfs:NAR.rgg\&quot;,\&quot;pfs:Untitled-1\&quot;&quot;"/>
          <de.grogra.pf.registry.Option name="panelId" type="java.lang.String" value="/ui/panels/texteditor"/>
          <de.grogra.pf.registry.Option name="panelTitle" type="java.lang.String" value="jEdit - NAR.rgg"/>
          <de.grogra.pf.registry.Option name="selected" type="java.lang.String" value="pfs:NAR.rgg"/>
         </de.grogra.pf.ui.registry.PanelFactory>
         <de.grogra.pf.registry.Link source="/ui/panels/attributeeditor"/>
        </de.grogra.pf.ui.registry.Tab>
        <de.grogra.pf.ui.registry.Tab selectedIndex="1">
         <de.grogra.pf.registry.Link source="/ui/panels/log"/>
         <de.grogra.pf.registry.Link source="/ui/panels/rgg/console"/>
        </de.grogra.pf.ui.registry.Tab>
       </de.grogra.pf.ui.registry.Split>
      </de.grogra.pf.ui.registry.Split>
     </de.grogra.pf.ui.registry.MainWindow>
     <de.grogra.pf.ui.registry.FloatingWindow height="475" width="696">
      <de.grogra.pf.ui.registry.PanelFactory source="/ui/panels/chart">
       <de.grogra.pf.registry.Option name="dataset" type="java.lang.String" value="/project/objects/datasets/dat"/>
       <de.grogra.pf.registry.Option name="panelId" type="java.lang.String" value="/ui/panels/chart?dat"/>
       <de.grogra.pf.registry.Option name="panelTitle" type="java.lang.String" value="dat"/>
       <de.grogra.pf.registry.Option name="plot" type="java.lang.Integer" value="7"/>
      </de.grogra.pf.ui.registry.PanelFactory>
     </de.grogra.pf.ui.registry.FloatingWindow>
    </de.grogra.pf.ui.registry.Layout>
   </ref>
  </ref>
 </registry>
</project>
