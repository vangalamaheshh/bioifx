Ext.onReady(buildTree);

function buildTree() {
	
	var treeLoader = new Ext.tree.TreeLoader({ 
		dataUrl: '../cgi-bin/NcbiJsonSender.cgi'
	});
	
	treeLoader.on("beforeload", function(treeLoader,node){
		treeLoader.baseParams.checked = node.attributes.checked;
	});
	
	var myStore = new Ext.data.JsonStore({
		url : '../cgi-bin/SendGridJson.cgi',
		root : 'info',
		fields : ['orgName', 'refseqId', {
			name : 'seqLen',
			type : 'int'
		}]
	});
	
	var clovrMug = new Ext.Button({
		text : 'CloVR-mugsy',
		handler : function() {
			var tree = Ext.getCmp('tree-container');
			var nodeIds = '',selNodes = tree.getChecked();
			Ext.each(selNodes, function(node){
				if(nodeIds.length > 0){
					nodeIds += ', ';
				}
				nodeIds += node.id;
			});
       	
			if(!nodeIds){
				Ext.MessageBox.show({
					title : 'Message',
					msg : 'You have not selected any ref seq, ' + 
						'please select the organims of interest in the tree widget beside',
					icon : Ext.MessageBox.INFO,
					buttons : Ext.Msg.OK,
					closable : false
		   		});
			}
			else {
				Ext.Ajax.request({
					url : '../cgi-bin/NcbiJsonSender.cgi',
					form : 'clovr_mug',
					params : {
						IDs : nodeIds,
						pipeline : 'clovr_mugsy'
					},
					success : function(response) {
						Ext.MessageBox.alert('Message', 'Pipeline invocation success. Task : ' + response.responseText);
					},
					failure : function(response) {
						Ext.MessageBox.alert('Message', 'Pipeline invocation failed');
					}
				});
			}
		}		
	});
	
	var clovrJoc = new Ext.Button({
		text : 'CloVR-JOC',
		handler : function() {
			var tree = Ext.getCmp('tree-container');
			var nodeIds = '',selNodes = tree.getChecked();
			Ext.each(selNodes, function(node){
				if(nodeIds.length > 0){
					nodeIds += ', ';
				}
				nodeIds += node.id;
			});
       		
			if(!nodeIds){
				Ext.MessageBox.show({
					title : 'Message',
					msg : 'You have not selected any ref seq, ' + 
						'please select the organims of interest in the tree widget beside',
					icon : Ext.MessageBox.INFO,
					buttons : Ext.Msg.OK,
					closable : false
	      		});
			} 
			else {
				Ext.Ajax.request({
					url : '../cgi-bin/NcbiJsonSender.cgi',
					form : 'clovr_joc',
					params : {
						IDs : nodeIds,
						pipeline : 'clovr_comparative'
					},
					success : function(response) {
						Ext.MessageBox.alert('Message', 'Pipeline invocation success. Task : ' + response.responseText);
					},
					failure : function(response) {
						Ext.MessageBox.alert('Message', 'Pipeline invocation failed');
					}
				});
			}
		}
	});
	
	var clovrPan = new Ext.Button({
		text : 'CloVR-pangenome',
		handler : function() {
			var tree = Ext.getCmp('tree-container');
			var nodeIds = '',selNodes = tree.getChecked();
			Ext.each(selNodes, function(node){
				if(nodeIds.length > 0){
					nodeIds += ', ';
				}
				nodeIds += node.id;
			});
   		
			if(!nodeIds){
				Ext.MessageBox.show({
					title : 'Message',
					msg : 'You have not selected any ref seq, ' + 
						'please select the organims of interest in the tree widget beside',
					icon : Ext.MessageBox.INFO,
					buttons : Ext.Msg.OK,
					closable : false
	       		});
			}
			
			else {
				Ext.Ajax.request({
					url : '../cgi-bin/NcbiJsonSender.cgi',
					form : 'clovr_pan',
					params : {
						IDs : nodeIds,
						pipeline : 'clovr_pangenome'
					},
					success : function(response) {
						Ext.MessageBox.alert('Message', 'Pipeline invocation success. Task : ' + response.responseText);
					},
					failure : function(response) {
						Ext.MessageBox.alert('Message', 'Pipeline invocation failed');
					}
				});
			}
		}
	});
	
	var myBar = new Ext.Toolbar({
		items : [{
			xtype : 'button',
			text : 'CloVR-comparative',
			handler : function() {
				Ext.Msg.alert('Message','clovr comparative is under construction');
			}
		}, '->', clovrPan, '-', clovrJoc, '-', clovrMug]
	});
	
	var myGrid = new Ext.grid.GridPanel({
		title : 'User Selected Ncbi Refseq Information',
		store : myStore,
		renderTo : Ext.get('mainPanel'),
		height : 500,
		//width : 800,
		columns : [ new Ext.grid.RowNumberer(), {
			header : 'ID',
			width : 30,
			dataIndex : 'id',
			sortable : true,
			hidden : true
		}, {
			id : 'org-name',
			header : 'Organism name', 
			width : 100,
			dataIndex : 'orgName',
			sortable : true
		}, {
			header : 'Sequence Length',
			width : 150,
			dataIndex : 'seqLen',
			sortable : true
		}, {
			header : 'Refseq ID',
			width : 150,
			dataIndex : 'refseqId',
			sortable : true,
			align : 'center'
		}],
		autoExpandColumn : 'org-name',
		loadMask : true,
		columnLines : true,
		bbar : myBar
	});
	
	var comboStore = new Ext.data.ArrayStore({
		url : '../cgi-bin/GetLineage.cgi',
		id : 0,
		fields : ['nodeName']
	});
	
	var combo = new Ext.form.ComboBox({
		xtype : 'combobox',
		store : comboStore,
		width : 300,
		emptyText : 'Type a node name to auto display',
		displayField : 'nodeName',
		valueField : 'nodeName',
		editable : true,
		lazyRender : true,
		mode : 'remote',
		forceSelection : true,
		triggerAction : 'all',
		typeAhead : true,
		typeAheadDelay : 300, 
	});
	
	comboStore.load();
	
	combo.on('select', function() {
		Ext.Ajax.request({
			url : '../cgi-bin/GetLineage.cgi',
			method : 'POST',
			params : {
				selectedNode : combo.getValue()
			},
			success : function(response) {
				var treeContainer = Ext.getCmp('tree-container');				
				var array = Ext.util.JSON.decode(response.responseText);
				var node;
				Ext.each(array, function(item, index, allItems) {
					node = treeContainer.getNodeById(item);
					node.expand();
				});
				node.select();
				//node.fireEvent('click');
				myStore.load({
					params : {
						id : node.attributes.id
					}
				});
			}
		});
	});
	
	new Ext.Viewport({
		layout : 'border',
		items : [{
			region : 'west',
		   	title : 'Navigation',
		   	id : 'tree-container',
		   	width : 300,
		   	xtype : 'treepanel',
		   	tbar : [combo],
		   	autoScroll : true,
		   	collapsible : true,
		   	split : true,
			useArrows : true,
		   	loader : treeLoader,
		   	root : new Ext.tree.AsyncTreeNode({
				text : 'Root Node',
				id : '/',
		   		checked : false
		   	}),
		   	rootVisible : false,
		   	listeners : {
				'click' : function(node) {
					myStore.load({
						params : {
							id : node.attributes.id
						}
					});
				},
				'checkchange' : function(node, checked){
					node.eachChild(function(n) {
		    			n.getUI().toggleCheck(checked);
		    		});
		    	}
		    } 
	 	}, {
	 			region : 'center',
	 			title : 'Clovr Comparative',
				layout : 'fit',
	 			id : 'mainPanel',
	 			items : [myGrid]
		}]
	});
}
