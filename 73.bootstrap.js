"use strict";(self.webpackChunk=self.webpackChunk||[]).push([[73],{83:(t,e,n)=>{n.d(e,{A:()=>a});var s=n(354),i=n.n(s),r=n(314),o=n.n(r)()(i());o.push([t.id,"/* TODO: Minimize CSS: https://webpack.js.org/plugins/mini-css-extract-plugin/#minimizing-for-production */\n\n* {\n    font-family: sans-serif;\n}\n\nhtml,\nbody {\n    margin: 0;\n    padding: 0;\n    overflow: hidden;\n}\n\nbody {\n    background: linear-gradient(0.25turn, #BBCCEE, #CCEEFF, #CCDDAA, #EEEEBB, #FFCCCC);\n}\n\nheader {\n    height: 2em;\n    width: 100%;\n    background: #FFF8;\n    display: flex;\n    gap: 1em;\n    justify-content: center;\n    align-items: center;\n}\n\n#svg-puzzle-container,\n#svg-puzzle {\n    height: 100%;\n    width: 100%;\n    max-width: 100vw;\n    max-height: calc(100vh - 2em);\n}\n\n#svg-puzzle>path {\n    stroke-width: 0.025;\n}\n\n#svg-puzzle>path.dragging {\n    filter: contrast(1.5);\n}\n\n#svg-puzzle>path.shadowGoal {\n    fill-opacity: 0.25;\n    stroke-opacity: 0.5;\n}\n\n.puzzle-string {\n    font-family: monospace;\n}","",{version:3,sources:["webpack://./style.css"],names:[],mappings:"AAAA,0GAA0G;;AAE1G;IACI,uBAAuB;AAC3B;;AAEA;;IAEI,SAAS;IACT,UAAU;IACV,gBAAgB;AACpB;;AAEA;IACI,kFAAkF;AACtF;;AAEA;IACI,WAAW;IACX,WAAW;IACX,iBAAiB;IACjB,aAAa;IACb,QAAQ;IACR,uBAAuB;IACvB,mBAAmB;AACvB;;AAEA;;IAEI,YAAY;IACZ,WAAW;IACX,gBAAgB;IAChB,6BAA6B;AACjC;;AAEA;IACI,mBAAmB;AACvB;;AAEA;IACI,qBAAqB;AACzB;;AAEA;IACI,kBAAkB;IAClB,mBAAmB;AACvB;;AAEA;IACI,sBAAsB;AAC1B",sourcesContent:["/* TODO: Minimize CSS: https://webpack.js.org/plugins/mini-css-extract-plugin/#minimizing-for-production */\r\n\r\n* {\r\n    font-family: sans-serif;\r\n}\r\n\r\nhtml,\r\nbody {\r\n    margin: 0;\r\n    padding: 0;\r\n    overflow: hidden;\r\n}\r\n\r\nbody {\r\n    background: linear-gradient(0.25turn, #BBCCEE, #CCEEFF, #CCDDAA, #EEEEBB, #FFCCCC);\r\n}\r\n\r\nheader {\r\n    height: 2em;\r\n    width: 100%;\r\n    background: #FFF8;\r\n    display: flex;\r\n    gap: 1em;\r\n    justify-content: center;\r\n    align-items: center;\r\n}\r\n\r\n#svg-puzzle-container,\r\n#svg-puzzle {\r\n    height: 100%;\r\n    width: 100%;\r\n    max-width: 100vw;\r\n    max-height: calc(100vh - 2em);\r\n}\r\n\r\n#svg-puzzle>path {\r\n    stroke-width: 0.025;\r\n}\r\n\r\n#svg-puzzle>path.dragging {\r\n    filter: contrast(1.5);\r\n}\r\n\r\n#svg-puzzle>path.shadowGoal {\r\n    fill-opacity: 0.25;\r\n    stroke-opacity: 0.5;\r\n}\r\n\r\n.puzzle-string {\r\n    font-family: monospace;\r\n}"],sourceRoot:""}]);const a=o},804:(t,e,n)=>{var s=n(72),i=n.n(s),r=n(825),o=n.n(r),a=n(659),l=n.n(a),c=n(56),_=n.n(c),h=n(540),u=n.n(h),g=n(113),f=n.n(g),d=n(83),m={};m.styleTagTransform=f(),m.setAttributes=_(),m.insert=l().bind(null,"head"),m.domAPI=o(),m.insertStyleElement=u(),i()(d.A,m),d.A&&d.A.locals&&d.A.locals},73:(t,e,n)=>{n.a(t,(async(t,s)=>{try{n.r(e),n(804);var i=n(410),r=n(617),o=n.n(r),a=n(883),l=t([i]);i=(l.then?(await l)():l)[0];const c=new a.A,_=16,h=(1<<_)-1,u=1<<_-1,g=32;function f(t,e){const n=(Math.round(t*g)&h)>>>0,s=(Math.round(e*g)&h)>>>0;return n<<_|s}function d(t){const e=t>>>_&h;return(e&u?e|~h:e)/g}function m(t){const e=t&h;return(e&u?e|~h:e)/g}const p=(t,e,n)=>f(e+d(t),n+m(t)),A=(t,e)=>f(d(t)+d(e),m(t)+m(e)),b=(t,e)=>A(t,w(e,-1)),w=(t,e)=>f(d(t)*e,m(t)*e),v=(t,e)=>new Set([...t].map((t=>A(t,e)))),k=(t,e)=>v(t,w(e,-1)),y=(t,e)=>d(t)*d(e)+m(t)*m(e),C=t=>t*t,z=(t,e)=>Math.sqrt(C(d(t)-d(e))+C(m(t)-m(e)));function E(t,e){return[...t].every((t=>!e.has(t)))}function B(t){let e=new Set;for(const n of t)for(const t of n)e.add(t);return e}function I(t){let e=Number.MAX_SAFE_INTEGER,n=Number.MAX_SAFE_INTEGER,s=Number.MIN_SAFE_INTEGER,i=Number.MIN_SAFE_INTEGER;for(const r of t)s=Math.max(s,d(r)),i=Math.max(i,m(r)),e=Math.min(e,d(r)),n=Math.min(n,m(r));return[f(e,n),f(s,i)]}function x(t,e){var n;const s=[f(1,0),f(0,-1),f(-1,0),f(0,1)];let i=new Map;t.forEach((e=>{const n=p(e,-1,0),s=p(e,1,0),r=p(e,0,-1),o=p(e,0,1);t.has(n)||i.set(p(e,-.5,0),1),t.has(s)||i.set(p(e,.5,0),3),t.has(r)||i.set(p(e,0,-.5),0),t.has(o)||i.set(p(e,0,.5),2)}));let r,o="",a=0;for(;(r=null===(n=i.keys().next())||void 0===n?void 0:n.value)&&a++<1e4;){let n=r,l=i.get(n),c=A(n,w(s[(l+3)%4],e));o+=`M${d(c)} ${m(c)}`;do{const r=s[l],a=w(r,.5),_=(w(r,-.5),w(s[(l+1)%4],.5)),h=w(s[(l+3)%4],.5),u=A(n,a),g=s[(l+1)%4],f=r;t.has(A(u,A(_,a)))?(n=A(u,_),l=(l+1)%4):t.has(A(u,A(h,a)))?n=A(u,a):(n=A(u,h),l=(l+3)%4);const p=c;c=A(n,w(s[(l+3)%4],e));const v=A(u,A(w(g,y(g,b(p,u))),w(f,y(f,b(c,u)))));o+=`C${d(v)} ${m(v)} ${d(v)} ${m(v)} ${d(c)} ${m(c)}`,i.delete(n)}while(n!=r&&a++<1e4);o+="Z"}return o}const S=".";function $(t,e=.75,n=1){return t==S?"#BBB":`oklch(${e} 0.2 ${65557*t.charCodeAt(0)%360} / ${n})`}const M="http://www.w3.org/2000/svg",F=document.getElementById("svg-puzzle-container"),j=document.createElementNS(M,"svg");j.id="svg-puzzle",j.setAttribute("xmlns",M),F.appendChild(j);class N{constructor(t,e){let[n,s]=I(t);this.shape=k(t,n),this.offset=n,this.char=e}client_to_point(t,e){this.svg_point.x=t,this.svg_point.y=e;const n=this.svg_point.matrixTransform(this.svg_sctm_inverse);return f(n.x,n.y)}update_walking_offsets(t){this.walking_offsets_time+=t.length,this.walking_offsets=t.concat(this.walking_offsets)}update_translation(t){const e=this.walking_offsets,n=this.walking_offsets_time;if(e&&e.length>0){const s=Math.max(0,n-t*Math.pow(Math.ceil(n),.75)/100);for(this.walking_offsets_time=s;e.length>s+2;)e.pop();if(1===e.length){const t=this.walking_offsets[0];this.path.style.transform=`translate(${d(t)}px, ${m(t)}px)`}else{const t=s%1,e=this.walking_offsets[Math.ceil(s)],n=this.walking_offsets[Math.floor(s)];this.path.style.transform=`translate(${t*d(e)+(1-t)*d(n)}px, ${t*m(e)+(1-t)*m(n)}px)`}}}get_coordinates(){return v(this.shape,this.offset)}construct_elem_path(t=1/32){let e=document.createElementNS(M,"path");if(e.setAttribute("d",x(this.shape,this.char===S?-t:t)),e.setAttribute("fill-rule","evenodd"),e.setAttribute("stroke",$(this.char,.25)),this.char!==S){const t=`block-pattern-${this.char}`;e.setAttribute("fill",`url(#${t})`)}else e.setAttribute("fill","#CCC");return e}create_elem_pattern(t){var e;if(this.char!==S){const n=null!==(e=t.querySelector("defs"))&&void 0!==e?e:t.appendChild(document.createElementNS(M,"defs")),s=document.createElementNS(M,"pattern"),i=`block-pattern-${this.char}`,[r,o]=I(this.shape);s.setAttribute("id",i),s.setAttribute("patternUnits","userSpaceOnUse"),s.setAttribute("x","-0.5"),s.setAttribute("y","-0.5"),s.setAttribute("width",`${d(o)+1}`),s.setAttribute("height",`${m(o)+1}`);const a=document.createElementNS(M,"rect");a.setAttribute("width",`${d(o)+1}`),a.setAttribute("height",`${m(o)+1}`),a.setAttribute("fill",$(this.char)),s.appendChild(a);const l=this.char.charCodeAt(0);for(let t of this.shape)for(let e=0;e<10;e++){const n=document.createElementNS(M,"text"),i=(52.09281985113131*l+13.087032255978896*e)%1+d(t),r=(28.640673508054988*l+94.82420753083849*e)%1+m(t),o=(14.33696587126313*l+40.125163576904164*e)%360,a=(72.31364428954059*l+61.41388485532069*e)%.25,c=(75.42795940481496*l+85.34675348929252*e)%.4+.2;n.setAttribute("x",i.toString()),n.setAttribute("y",r.toString()),n.setAttribute("fill",$(this.char,.5,a)),n.setAttribute("font-size",c.toString()),n.setAttribute("transform",`rotate(${o} ${i} ${r})`),n.textContent=this.char,s.appendChild(n)}n.appendChild(s)}}initialise_elem(t){this.create_elem_pattern(t);const e=this.construct_elem_path();this.path=e,this.svg_elem=t,this.svg_sctm_inverse=t.getScreenCTM().inverse(),this.svg_point=t.createSVGPoint(),t.appendChild(this.path),this.walking_offsets=[],this.walking_offsets_time=0,this.update_walking_offsets([this.offset]),this.update_translation(0)}make_interactive(t,e){if(!this.path||!this.svg_elem)return;const n=t.blockstate,s=n.bounds.get_coordinates();let i,r,a=this.offset;o()(this.path).draggable({listeners:{start:t=>{a=this.offset,i=b(this.client_to_point(t.client.x,t.client.y),a),this.path.classList.add("dragging");const o=B(n.blocks.filter(((t,n)=>n!=e)).map((t=>t.get_coordinates()))),l=t=>{const e=v(this.shape,t);return n=e,i=s,!![...n].every((t=>i.has(t)))&&!!E(e,o);var n,i};let c=[a];for(r=new Set([a]);c.length>0;){let t=c.shift();for(let e of[A(t,f(0,1)),A(t,f(0,-1)),A(t,f(1,0)),A(t,f(-1,0))])!r.has(e)&&l(e)&&(r.add(e),c.push(e))}},move:t=>{const e=this.offset,n=b(this.client_to_point(t.client.x,t.client.y),i);let s=a,o=1/0;for(let t of r){const e=z(n,t);e<o&&(o=e,s=t)}if(e!=s){let t=[e],n={old_offset:null};for(;t.length>0;){let e=t.shift();if(e==s)break;for(let s of[A(e,f(0,1)),A(e,f(0,-1)),A(e,f(1,0)),A(e,f(-1,0))])!(s in n)&&r.has(s)&&(n[s]=e,t.push(s))}let i=[s],o=n[s];for(;o!==e;)i.push(o),o=n[o];this.offset=s,this.update_walking_offsets(i)}},end:e=>{if(a!==this.offset&&(t.move_counter+=1,U.textContent=t.move_counter.toString(),t.won())){const e=24;null!==t.min_moves&&t.move_counter<t.min_moves?c.addConfetti({emojis:["🐞"],confettiNumber:e}):null!==t.min_moves&&t.move_counter===t.min_moves?c.addConfetti({emojis:["🏆"],confettiNumber:e}):c.addConfetti({confettiNumber:e,confettiColors:t.blockstate.blocks.map((t=>$(t.char)))})}this.path.classList.remove("dragging")}}})}}class T{constructor(t,e){[this.min,this.max]=I(t.shape),this.bounds=t,this.blocks=e}static blockstate_from_string(t){const e=t.split(/\r?\n/g);let n=Number.MAX_SAFE_INTEGER,s=Number.MAX_SAFE_INTEGER,i={},r=new Set,o=0;for(let t of e){let e=0;for(let a of[...t])/\s/.test(a)||(n=Math.min(n,e),s=Math.min(s,o),a!=S&&(i[a]=i[a]||new Set,i[a].add(f(e,o))),r.add(f(e,o))),e++;o++}const a=f(n,s);let l=new N(k(r,a),S),c=Object.entries(i).map((([t,e])=>new N(k(e,a),t)));return new T(l,c)}to_string(){let t="";const e=this.blocks.map((t=>[t,t.get_coordinates()])),n=this.bounds.get_coordinates();for(let s=m(this.min);s<=m(this.max);s++){for(let i=d(this.min);i<=d(this.max);i++){let r=!1;for(let[n,o]of e)if(o.has(f(i,s))){t+=n.char,r=!0;break}r||(n.has(f(i,s))?t+=S:t+=" ")}t+="\n"}return t}initialise(t,e){t.innerHTML="";const n=d(this.max)-d(this.min)+1,s=m(this.max)-m(this.min)+1;t.setAttribute("viewBox",`${d(this.min)-1} ${m(this.min)-1} ${n+1} ${s+1}`),this.svg_elem=t,this.bounds.initialise_elem(t);for(let n=0;n<this.blocks.length;n++)if(null!==e[n]){const s=this.blocks[n].construct_elem_path(0),i=e[n];s.classList.add("shadowGoal"),s.setAttribute("transform",`translate(${d(i)}, ${m(i)})`),t.appendChild(s)}for(let e of this.blocks)e.initialise_elem(t)}make_interactive(t){for(let e=0;e<this.blocks.length;e++)this.blocks[e].make_interactive(t,e)}}class G{constructor(t,e,n=null){this.min_moves=n,this.move_counter=0;const s=T.blockstate_from_string(t);this.blockstate=s;const i=T.blockstate_from_string(e);this.goal_offsets=this.blockstate.blocks.map((t=>{for(let e of i.blocks)if(e.char===t.char)return e.offset;return null})),this.start_string=s.to_string(),this.goal_string=i.to_string()}won(){const t=this.goal_offsets;return this.blockstate.blocks.every(((e,n)=>null===t[n]||e.offset===t[n]))}initialise(t){this.blockstate=T.blockstate_from_string(this.start_string),this.blockstate.initialise(t,this.goal_offsets),this.blockstate.make_interactive(this),this.move_counter=0,U.textContent="0"}}const L=document.getElementById("change-puzzle-dialog"),D=document.getElementById("change-puzzle-btn"),R=document.getElementById("reset-puzzle-btn"),q=document.getElementById("puzzle-textarea-start"),O=document.getElementById("puzzle-textarea-goal"),X=document.getElementById("puzzle-submit-btn"),U=document.getElementById("move-counter");D.addEventListener("click",(()=>{O.value=P.goal_string,q.value=P.start_string,L.showModal()})),R.addEventListener("click",(()=>{P.initialise(j)})),X.addEventListener("click",(t=>{t.preventDefault();const e=q.value,n=O.value;L.close(),P=new G(e,n),P.initialise(j)}));const W=document.getElementById("puzzle-selection"),J=(0,i.J9)();for(let Y of J){const Z=document.createElement("option");Z.textContent=Y.name,Z.value=Y.name,W.appendChild(Z)}W.addEventListener("change",(()=>{const t=W.value,e=J.find((e=>e.name===t));q.value=T.blockstate_from_string(e.start).to_string(),O.value=T.blockstate_from_string(e.goal).to_string()}));let Q,P=new G(J[0].start,J[0].goal);function V(t){const e=Q?t-Q:0;Q=t;for(let t of P.blockstate.blocks)t.update_translation(e);requestAnimationFrame(V)}P.initialise(j),requestAnimationFrame(V),s()}catch(H){s(H)}}))},410:(t,e,n)=>{n.a(t,(async(t,s)=>{try{n.d(e,{J9:()=>r.J9});var i=n(631),r=n(518),o=t([i]);i=(o.then?(await o)():o)[0],(0,r.lI)(i),s()}catch(t){s(t)}}))},518:(t,e,n)=>{let s;function i(t){s=t}n.d(e,{J9:()=>w,Mk:()=>y,Qn:()=>z,lI:()=>i,yc:()=>C});let r=new("undefined"==typeof TextDecoder?(0,module.require)("util").TextDecoder:TextDecoder)("utf-8",{ignoreBOM:!0,fatal:!0});r.decode();let o=null;function a(){return null!==o&&0!==o.byteLength||(o=new Uint8Array(s.memory.buffer)),o}function l(t,e){return t>>>=0,r.decode(a().subarray(t,t+e))}const c=new Array(128).fill(void 0);c.push(void 0,null,!0,!1);let _=c.length;function h(t){_===c.length&&c.push(c.length+1);const e=_;return _=c[e],c[e]=t,e}let u=0,g=new("undefined"==typeof TextEncoder?(0,module.require)("util").TextEncoder:TextEncoder)("utf-8");const f="function"==typeof g.encodeInto?function(t,e){return g.encodeInto(t,e)}:function(t,e){const n=g.encode(t);return e.set(n),{read:t.length,written:n.length}};function d(t,e,n){if(void 0===n){const n=g.encode(t),s=e(n.length,1)>>>0;return a().subarray(s,s+n.length).set(n),u=n.length,s}let s=t.length,i=e(s,1)>>>0;const r=a();let o=0;for(;o<s;o++){const e=t.charCodeAt(o);if(e>127)break;r[i+o]=e}if(o!==s){0!==o&&(t=t.slice(o)),i=n(i,s,s=o+3*t.length,1)>>>0;const e=a().subarray(i+o,i+s);o+=f(t,e).written,i=n(i,s,o,1)>>>0}return u=o,i}let m=null;function p(){return null!==m&&0!==m.byteLength||(m=new Int32Array(s.memory.buffer)),m}let A=null;function b(t){const e=function(t){return c[t]}(t);return function(t){t<132||(c[t]=_,_=t)}(t),e}function w(){try{const i=s.__wbindgen_add_to_stack_pointer(-16);s.get_all_js_examples(i);var t=p()[i/4+0],e=p()[i/4+1],n=function(t,e){t>>>=0;const n=(null!==A&&0!==A.byteLength||(A=new Uint32Array(s.memory.buffer)),A).subarray(t/4,t/4+e),i=[];for(let t=0;t<n.length;t++)i.push(b(n[t]));return i}(t,e).slice();return s.__wbindgen_free(t,4*e,4),n}finally{s.__wbindgen_add_to_stack_pointer(16)}}const v="undefined"==typeof FinalizationRegistry?{register:()=>{},unregister:()=>{}}:new FinalizationRegistry((t=>s.__wbg_jspuzzle_free(t>>>0)));class k{static __wrap(t){t>>>=0;const e=Object.create(k.prototype);return e.__wbg_ptr=t,v.register(e,e.__wbg_ptr,e),e}__destroy_into_raw(){const t=this.__wbg_ptr;return this.__wbg_ptr=0,v.unregister(this),t}free(){const t=this.__destroy_into_raw();s.__wbg_jspuzzle_free(t)}get name(){let t,e;try{const r=s.__wbindgen_add_to_stack_pointer(-16);s.__wbg_get_jspuzzle_name(r,this.__wbg_ptr);var n=p()[r/4+0],i=p()[r/4+1];return t=n,e=i,l(n,i)}finally{s.__wbindgen_add_to_stack_pointer(16),s.__wbindgen_free(t,e,1)}}set name(t){const e=d(t,s.__wbindgen_malloc,s.__wbindgen_realloc),n=u;s.__wbg_set_jspuzzle_name(this.__wbg_ptr,e,n)}get start(){let t,e;try{const r=s.__wbindgen_add_to_stack_pointer(-16);s.__wbg_get_jspuzzle_start(r,this.__wbg_ptr);var n=p()[r/4+0],i=p()[r/4+1];return t=n,e=i,l(n,i)}finally{s.__wbindgen_add_to_stack_pointer(16),s.__wbindgen_free(t,e,1)}}set start(t){const e=d(t,s.__wbindgen_malloc,s.__wbindgen_realloc),n=u;s.__wbg_set_jspuzzle_start(this.__wbg_ptr,e,n)}get goal(){let t,e;try{const r=s.__wbindgen_add_to_stack_pointer(-16);s.__wbg_get_jspuzzle_goal(r,this.__wbg_ptr);var n=p()[r/4+0],i=p()[r/4+1];return t=n,e=i,l(n,i)}finally{s.__wbindgen_add_to_stack_pointer(16),s.__wbindgen_free(t,e,1)}}set goal(t){const e=d(t,s.__wbindgen_malloc,s.__wbindgen_realloc),n=u;s.__wbg_set_jspuzzle_goal(this.__wbg_ptr,e,n)}get min_moves(){return s.__wbg_get_jspuzzle_min_moves(this.__wbg_ptr)>>>0}set min_moves(t){s.__wbg_set_jspuzzle_min_moves(this.__wbg_ptr,t)}}function y(t){return h(k.__wrap(t))}function C(t,e){return h(l(t,e))}function z(t,e){throw new Error(l(t,e))}},631:(t,e,n)=>{var s=n(518);t.exports=n.v(e,t.id,"3920994c0ac15281eb8a",{"./sliding_blocks_bg.js":{__wbg_jspuzzle_new:s.Mk,__wbindgen_string_new:s.yc,__wbindgen_throw:s.Qn}})}}]);
//# sourceMappingURL=73.bootstrap.js.map